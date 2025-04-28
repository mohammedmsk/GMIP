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
data=combined_results[,c("GWAS_name", "feature_name", "Method_name", "CV_strategy", "adj_p_value", "tau_normalized_original")]

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Function to generate and save individual heatmaps for each GWAS and CV strategy
generate_and_save_heatmap=function(gwas, strategy) {
  
  # Subset the data
  subset_data=data %>% 
    filter(GWAS_name == gwas, CV_strategy == strategy) %>%
    select(feature_name, Method_name, adj_p_value, tau_normalized_original)
  
  # Create matrices for p-values and tau-normalized values
  tau_matrix=reshape2::acast(subset_data, feature_name ~ Method_name, value.var="tau_normalized_original")
  p_matrix=reshape2::acast(subset_data, feature_name ~ Method_name, value.var="adj_p_value")
  
  # Log-transform p-values for better visualization
  p_matrix_log=-log10(p_matrix)
  
  # Define heatmap color scales with reversed colors for p-values
  tau_range=range(tau_matrix, na.rm=TRUE)
  col_fun_tau=colorRamp2(seq(tau_range[1], tau_range[2], length.out=3), c("blue", "white", "red"))
  
  p_range=range(p_matrix_log, na.rm=TRUE)
  col_fun_p=colorRamp2(seq(p_range[1], p_range[2], length.out=3), c("black", "white", "gold"))
  
  # Create the heatmaps without clustering
  heatmap_tau=Heatmap(tau_matrix, 
                         name="Tau Normalized", 
                         col=col_fun_tau, 
                         row_title="Features", 
                         column_title="Methods", 
                         heatmap_legend_param=list(title="Tau Normalized"),
                         na_col="grey",
                         cluster_rows=TRUE,   # Disable row clustering
                         cluster_columns=FALSE # Disable column clustering
                         )
  
  heatmap_p=Heatmap(p_matrix_log, 
                       name="-log10(p-value)", 
                       col=col_fun_p, 
                       row_title="Features", 
                       column_title="Methods", 
                       heatmap_legend_param=list(title="-log10(p-value)"),
                       na_col="grey",
                       cluster_rows=TRUE,    # Disable row clustering
                       cluster_columns=FALSE  # Disable column clustering
                       )
  
  # Combine the heatmaps side by side
  combined_heatmap=heatmap_tau + heatmap_p
  
  # Save the individual heatmap to a high-resolution PNG file
  png(filename=paste0("heatmaps/", gwas, "_", strategy, "_heatmap.png"), 
      width=2400, height=1600, res=300) # Adjust resolution for high-quality output
  draw(combined_heatmap, column_title=paste("GWAS:", gwas, "| CV Strategy:", strategy), 
       column_title_gp=gpar(fontsize=14, fontface="bold"))
  dev.off()
}

# Ensure the "heatmaps" directory exists
if (!dir.exists("heatmaps")) {
  dir.create("heatmaps")
}

# Generate and save heatmaps for each GWAS and CV strategy combination as high-resolution PNGs
unique_gwas=unique(data$GWAS_name)
unique_strategies=unique(data$CV_strategy)

for (gwas in unique_gwas) {
  for (strategy in unique_strategies) {
    generate_and_save_heatmap(gwas, strategy)
  }
}



