# Load necessary libraries
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer) # For ColorBrewer palettes
library(patchwork)    # For combining ggplot objects

base_dir="/scratch/04179/mshabb/GMIP2/2024_09_27"
subdirs=c(
  "NAFLD/expression_pops",  "NAFLD/full_pops", "NAFLD/nafld", "NAFLD/naga", "NAFLD/netwas", "NAFLD/pathway_pops", "NAFLD/ppi_pops",
  "RA/expression_pops",  "RA/full_pops", "RA/nafld", "RA/naga", "RA/netwas", "RA/pathway_pops", "RA/ppi_pops",
  "SCZ/expression_pops",  "SCZ/full_pops", "SCZ/nafld", "SCZ/naga", "SCZ/netwas", "SCZ/pathway_pops", "SCZ/ppi_pops")

mid_dirs="outdir/6_eval_results"
inner_dirs=c("naga_kfold", "naga_loco", "naga_nocv", "pops_kfold", "pops_loco", "pops_nocv") 
file_pattern="*.ypred_asRL.gsea.tsv"

# Collect all file paths and associated information
files_info <- rbindlist(lapply(subdirs, function(subdir) {
  file.list <- lapply(inner_dirs, function(inner_dir) {
    file_names <- list.files(
      path = file.path(base_dir, subdir, mid_dirs, inner_dir), 
      pattern = file_pattern, 
      full.names = TRUE
    )
    if (length(file_names) > 0) {
      data.table(subdir = subdir, inner_dir = inner_dir, file_path = file_names)
    }
  })
  rbindlist(file.list, fill = TRUE)
}))

# Read all files and add directory columns
combined_results <- rbindlist(lapply(files_info$file_path, function(file) {
  dt <- fread(file)
  dt[, `:=`(subdir = files_info$subdir[files_info$file_path == file],
            inner_dir = files_info$inner_dir[files_info$file_path == file])]
}), fill = TRUE)

# Extract GWAS name, feature name, Method name, and CV strategy
combined_results[, `:=`(
  GWAS_name = tstrsplit(subdir, "/", fixed = TRUE)[[1]],
  feature_name = tstrsplit(subdir, "/", fixed = TRUE)[[2]],
  Method_name = tstrsplit(inner_dir, "_", fixed = TRUE)[[1]],
  CV_strategy = tstrsplit(inner_dir, "_", fixed = TRUE)[[2]]
)]

# Update feature names and select relevant columns for visualization
combined_results$feature_name <- paste0(combined_results$feature_name, "_fset")
data=combined_results[,c("GWAS_name", "feature_name", "Method_name", "CV_strategy", "NES", "p.adjust", "ID")]
data=data[data$ID=="top_0500_gwas_genes"]

# Create directory for barplots if it doesn't exist
if (!dir.exists("barplots_gsea_top0500")) {
  dir.create("barplots_gsea_top0500")
}

# Generate a unified color palette for all unique feature names across all data
unique_features <- unique(data$feature_name)
num_features <- length(unique_features)
palette <- brewer.pal(min(max(num_features, 3), 9), "Set1")  # Adjust the color palette to the number of features
feature_colors <- setNames(palette, unique_features)  # Create a named vector for consistent mapping

# Function to generate and save combined bar plots for each GWAS with all CV strategies combined in one plot
generate_combined_barplot_per_gwas <- function(gwas) {
  
  # Get unique CV strategies for the current GWAS
  unique_strategies <- unique(data$CV_strategy)
  
  # Create an empty list to store plots for each CV strategy
  strategy_plots <- list()
  
  # Loop through each CV strategy to generate individual combined bar plots
  for (strategy in unique_strategies) {
    
    # Subset the data for the given GWAS and CV strategy
    subset_data <- data %>% 
      filter(GWAS_name == gwas, CV_strategy == strategy)
    
    # Get unique method names
    method_names <- unique(subset_data$Method_name)
    
    # Determine the maximum y-axis limit across all methods for consistent scaling
    y_max <- max(subset_data$NES, na.rm = TRUE) * 1.1  # Add 10% cushion to ensure visibility of asterisks
    
    # Create individual plots for each method within this CV strategy
    plots <- list()
    
    for (method in method_names) {
      method_data <- subset_data %>% filter(Method_name == method)
      
      p <- ggplot(method_data, aes(x = feature_name, y = NES)) +
        geom_segment(aes(xend = feature_name, yend = 0, color = feature_name), size = 15) +
        geom_hline(yintercept = 1.5, linetype = "dotted", color = "red", size = 1) +
        geom_text(aes(label = ifelse(p.adjust < 0.05, "*", "")), 
                  vjust = -0.5, size = 6, color = "black") +  # Add asterisks for significant values
        scale_color_manual(values = feature_colors) +  # Apply the consistent color palette
        labs(title = paste(method, "-", strategy), x = "", y = "gsea NES") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18, vjust = 0.5),
              axis.text.y = element_text(size = 18),
              plot.title = element_text(size = 24, face = "bold"),
              axis.title.y = element_text(size = 18),
              legend.position = "none") +  # Removing legend since it's repeated
        ylim(0, y_max)  # Apply consistent y-axis limits
      plots[[method]] <- p
    }
    
    # Combine all method plots for this strategy using patchwork
    combined_strategy_plot <- wrap_plots(plots, ncol = length(method_names))
    
    # Add the combined plot for this CV strategy to the list
    strategy_plots[[strategy]] <- combined_strategy_plot
  }
  
  # Combine all CV strategy plots into a single plot with 3 rows
  final_combined_plot <- wrap_plots(strategy_plots, ncol = 1)
  
  # Save the combined plot as a high-resolution PNG file
  ggsave(filename = paste0("barplots_gsea_top0500/", gwas, "_combined_barplot_all_strategies.png"), 
         plot = final_combined_plot, width = 5 * length(method_names), height = 18, dpi = 300)  # Adjust the height to accommodate all rows
}

# Generate and save combined bar plots for each GWAS with all CV strategies in one plot
unique_gwas <- unique(data$GWAS_name)

for (gwas in unique_gwas) {
  generate_combined_barplot_per_gwas(gwas)
}
