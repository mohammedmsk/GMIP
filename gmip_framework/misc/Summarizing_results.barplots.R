library(data.table)
library(stringr)

# Define the fixed parts of your directory structure
base_dir="/scratch/04179/mshabb/GMIP2/2024_09_26"
subdirs=c("RA/naga", "RA/netwas", "SCZ/naga", "SCZ/netwas")
subdirs=c("RA/expression_pops",  "RA/full_pops", "RA/nafld", "RA/naga", "RA/netwas", "RA/pathway_pops", "RA/ppi_pops",
  "SCZ/expression_pops",  "SCZ/full_pops", "SCZ/nafld", "SCZ/naga", "SCZ/netwas", "SCZ/pathway_pops", "SCZ/ppi_pops")
mid_dirs="outdir/5_bm_results"
inner_dirs=c("naga_kfold", "naga_loco", "naga_nocv", "pops_kfold", "pops_loco", "pops_nocv") 
final_dir="part_2"
file_pattern="*.results_withStats"

# Collect all file paths and associated information
files_info=rbindlist(lapply(subdirs, function(subdir) {
  file.list=lapply(inner_dirs, function(inner_dir) {
    file_names=list.files(
      path=file.path(base_dir, subdir, mid_dirs, inner_dir, final_dir), 
      pattern=file_pattern, 
      full.names=TRUE
    )
    if (length(file_names) > 0) {
      data.table(subdir=subdir, inner_dir=inner_dir, file_path=file_names)
    }
  })
  rbindlist(file.list, fill=TRUE)
}))

# Read all files and add directory columns
combined_results=rbindlist(lapply(files_info$file_path, function(file) {
  dt=fread(file)
  dt[, `:=`(subdir=files_info$subdir[files_info$file_path == file],
            inner_dir=files_info$inner_dir[files_info$file_path == file])]
}), fill=TRUE)

# Extract the required information using vectorized operations
combined_results[, `:=`(
  GWAS_name=tstrsplit(subdir, "/", fixed=TRUE)[[1]],
  feature_name=tstrsplit(subdir, "/", fixed=TRUE)[[2]],
  Method_name=tstrsplit(inner_dir, "_", fixed=TRUE)[[1]],
  CV_strategy=tstrsplit(inner_dir, "_", fixed=TRUE)[[2]]
)]


# Adjust p-values using Bonferroni correction and create a copy of the original tau values
combined_results[, `:=`(
  adj_p_value=p.adjust(Enrichment_p, method="bonferroni"),
  tau_normalized_original=tau_normalized
)]

# Set tau_normalized to NA where adjusted p-value is greater than or equal to 0.05
combined_results[adj_p_value >= 0.05, tau_normalized := NA]

# Summarize the results to find the best feature_name, Method_name, and CV_strategy for each GWAS based on max tau_normalized
best_results=combined_results[, .SD[which.max(tau_normalized)],  by=GWAS_name, 
  .SDcols=c("feature_name", "Method_name", "CV_strategy", "tau_normalized", "tau_normalized_original")
]

combined_results$feature_name=paste0(combined_results$feature_name, "_features")
data=combined_results[,c("GWAS_name", "feature_name", "Method_name", "CV_strategy", "adj_p_value", "tau_normalized_original", "tau_standard_error_normalized")]

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer) # For ColorBrewer palettes
library(patchwork)    # For combining ggplot objects

# Ensure the "barplots" directory exists
if (!dir.exists("barplots")) {
  dir.create("barplots")
}

# Function to generate and save bar plots for each GWAS and CV strategy combination
generate_combined_barplot <- function(gwas, strategy) {
  
  # Subset the data for the given combination of GWAS and CV strategy
  subset_data <- data %>% 
    filter(GWAS_name == gwas, CV_strategy == strategy)
  
  # Get unique method names
  method_names <- unique(subset_data$Method_name)
  
  # Determine the maximum y-axis limit across all methods for consistent scaling
  y_max <- max(subset_data$tau_normalized_original, na.rm = TRUE)
  
  # Generate a color palette based on the number of unique features
  num_features <- length(unique(subset_data$feature_name))
  palette <- brewer.pal(min(max(num_features, 3), 9), "Set1")  # Using Set1 palette with a limit from 3 to 9 colors
  
  # Create individual plots for each method
  plots <- list()
  
  for (method in method_names) {
    method_data <- subset_data %>% filter(Method_name == method)
    
    p <- ggplot(method_data, aes(x = feature_name, y = tau_normalized_original)) +
      geom_segment(aes(xend = feature_name, yend = 0, color = feature_name), size = 20) +
      scale_color_manual(values = palette) +
      labs(title = method, x = "Feature", y = "Tau Normalized") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "none") +  # Removing legend since it's repeated
      ylim(0, y_max)  # Apply consistent y-axis limits
     
    plots[[method]] <- p
  }
  
  # Combine all method plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = length(method_names))
  
  # Save the combined plot as a high-resolution PNG file
  ggsave(filename = paste0("barplots/", gwas, "_", strategy, "_combined_barplot.png"), 
         plot = combined_plot, width = 5 * length(method_names), height = 6, dpi = 300)
}

# Generate and save combined bar plots for each GWAS and CV strategy combination
unique_gwas <- unique(data$GWAS_name)
unique_strategies <- unique(data$CV_strategy)

for (gwas in unique_gwas) {
  for (strategy in unique_strategies) {
    generate_combined_barplot(gwas, strategy)
  }
}

