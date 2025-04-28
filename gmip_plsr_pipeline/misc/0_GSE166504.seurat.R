#conda activate seurat
# Load libraries
library(tidyverse)
library(data.table)
library(Seurat)
library(irlba)
library(Matrix)
library(future)
library(ggrastr)
library(tidytext)
library(matrixTests)
library(BuenColors)
source("~/GMIP/misc/utils.R")
setwd("/scratch/04179/mshabb/NAFLD_db/GSE166504")

metadata=as.data.frame(data.table::fread('GSE166504_cell_metadata.20220204.tsv.gz'))
rownames(metadata)=paste0(metadata$FileName, "_", metadata$CellID)

raw_counts=data.table::fread('GSE166504_cell_raw_counts.20220204.txt.gz')
mat=as.matrix(raw_counts, rownames=1)
#source("~/utils.R")

# Parameters
name="GSE166504"
number_pcs=30
vargenes=2000
clus_res=0.2

# Setup
dir.create(paste0("plots/", name), recursive=T)
dir.create(paste0("features/", name), recursive=T)

mat2=ConvertToENSGAndProcessMatrix(mat, 'mouse_symbol')
dim(mat2)
new_seurat_object=CreateSeuratObject(counts=mat2)
metadata=metadata[colnames(new_seurat_object),]
new_seurat_object@meta.data=metadata
so=new_seurat_object
so=NormalizeData(so, normalization.method="LogNormalize", scale.factor=1000000)
so=ScaleData(so, min.cells.to.block=1, block.size=500)
dim(so)

# Identify variable genes
so=FindVariableFeatures(so, nfeatures=vargenes)
so=RunPCA(so, npcs=100)
so=ProjectDim(so, do.center=T)
so=RunICA(so, nics=number_pcs)
so=ProjectDim(so, reduction="ica", do.center=T)
so=FindNeighbors(so, dims=1:number_pcs, nn.eps=0)
so=FindClusters(so, resolution=clus_res, n.start=100)
so=RunUMAP(so, dims=1:number_pcs, min.dist=0.4, n.epochs=500, n.neighbors=10, learning.rate=0.1, spread=2)

SaveGlobalFeatures(so, name)

Idents(object=so) <- "CellType"
clus <- unique(so@meta.data$CellType)
demarkers_pre_def <- WithinClusterFeatures(so, "CellType", clus, name, suffix = "_pre_def1")

Idents(object=so) <- "grp"
clus <- unique(so@meta.data$grp)
demarkers_pre_def <- WithinClusterFeatures(so, "grp", clus, name, suffix = "_pre_def2")

Idents(object=so)="seurat_clusters"
clus=levels(so@meta.data$seurat_clusters)
demarkers=WithinClusterFeatures(so, "seurat_clusters", clus, name)

keep=read.table(
"/scratch/04179/mshabb/05577fee5848687884crrerfefe845/refdir//static_gits/pops/example/data/utils/gene_annot_jun10.txt",
  sep="\t", header=T, stringsAsFactors=F, col.names=c("ENSG", "symbol", "chr", "start", "end", "TSS"))


########################################################################################################################
#conda activate seurat
library(data.table)
setwd("/scratch/04179/mshabb/NAFLD_db/GSE166504/features/GSE166504")

# Step 1: Read all files
files=list.files(pattern = "\\.txt.gz$", full.names = TRUE)
dt_list=lapply(files, fread)
names(dt_list)=files

# Step 2: Append filename to column names and add prefix "GSE166504"
dt_list=lapply(seq_along(dt_list), function(i) {
  dt=dt_list[[i]]
  # Extracting filename without extension
  file_name=gsub(".txt.gz", "", basename(files[i]))
  # Modify column names except for the first one
  setnames(dt, old = names(dt)[-1], new = paste("GSE166504", file_name, names(dt)[-1], sep = "__"))
  return(dt)
})
names(dt_list)=files

# Step 3: Merge all data.tables on the first column
merged_dt=Reduce(function(x, y) merge(x, y, by="ENSG", all=TRUE), dt_list)
dim(merged_dt)
cols_to_standardize=names(merged_dt)[-1]

# Step 4: Standardize the columns except the first one using scale()
merged_dt[, (cols_to_standardize) := lapply(.SD, scale), .SDcols=cols_to_standardize]

# Add a 'control' column where zero rows have 0 and others have 1
# Identify rows where all columns to be standardized are zero
zero_rows=apply(dt_list[["./average_expression.txt.gz"]][,-1], 1, function(x) sum(x)==0)
table(zero_rows)
control_dt=data.table(ENSG=dt_list[["./average_expression.txt.gz"]][[1]], GSE166504_control=as.numeric(zero_rows))
dim(merged_dt)
merged_dt=merge(merged_dt, control_dt, by="ENSG", all=TRUE)
dim(merged_dt)
merged_dt[[629]]
# Calculate mean and standard deviation for each column
column_means=sapply(merged_dt[, -1], mean, na.rm = TRUE)
column_sds=sapply(merged_dt[, -1], sd, na.rm = TRUE)
dim(merged_dt)

# Step 5: Output the final data.table
fwrite(merged_dt, "../../GSE166504.standardized.final_output.txt", sep = "\t")
########################################################################################################################
#Steps for munge_features etc
system("mkdir -p /scratch/04179/mshabb/05577fee5848687884crrerfefe845/refdir/2_pops_files/nafld/GSE166504/")
system("cp ../../GSE166504.standardized.final_output.txt /scratch/04179/mshabb/05577fee5848687884crrerfefe845/refdir/2_pops_files/nafld/GSE166504/")
setwd("/scratch/04179/mshabb/05577fee5848687884crrerfefe845/refdir/2_pops_files/nafld/GSE166504/")
dir.create("feature_in", recursive=T, showWarnings=F)
system("mv GSE166504.standardized.final_output.txt feature_in")
system("rm -rf GSE166504_gene_features_munged")
system("mkdir -p GSE166504_gene_features_munged")
system("/work2/04179/mshabb/ls6/miniforge3/envs/GMIP/bin/python ~/GMIP/bin/pops/munge_feature_directory.py --gene_annot_path /scratch/04179/mshabb/05577fee5848687884crrerfefe845/refdir/static_gits/pops/example/data/utils/gene_annot_jun10.txt --feature_dir feature_in --save_prefix GSE166504_gene_features_munged/pops_features --max_cols 100")
system("echo GSE166504_control > control.features")

########################################################################################################################
#load("R_image.rds")
Idents(object=so) <- "CellType"
p = DimPlot(so)
p = p + theme(
  legend.text = element_text(size = 24),          # Adjust legend text size
  legend.title = element_text(size = 24),         # Adjust legend title size
  legend.key.size = unit(1.5, "lines")            # Adjust legend key size
)

ggsave(plot=p, "~/umap_celltype.png", width=18, height=12)

Idents(object=so) <- "grp"
p = DimPlot(so)
ggsave(plot=p, "~/umap_grp.png", width=15, height=12)

topdegenes <- demarkers_pre_def %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

topdegenes <- topdegenes %>%  group_by(cluster) %>% slice(c(1,2), with_ties = F)
topdegenes.df <- bind_cols(data.frame(t(so@assays$RNA@scale.data[topdegenes$ENSG,])),data.frame(so@reductions$umap@cell.embeddings)) %>%  as.tibble()

p <- topdegenes.df %>%
  reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
  merge(., keep, by.x = "variable", by.y = "ENSG") %>%
  as.tibble() %>%
  dplyr::mutate(value = case_when(value > 3 ~ 3,
                                  value < -3 ~ -3,
                                  T ~ value)) %>%
  dplyr::rename("exprs" = value) %>%
  ggplot(., aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(aes(color = exprs), size = 1, raster.dpi = 100) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot() +
  facet_wrap(~symbol, ncol = 4)
ggsave(plot=p, "~/umap_degenes.png", width=15, height=12)
