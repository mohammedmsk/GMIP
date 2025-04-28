This README explains how the pops features were setup.

```{r}
setwd("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir")
pops_features=data.table::fread("2_pops_files/data_from_dropbox/full_pops/feature_in/PoPS.features.txt.gz")
dim(pops_features)
pops_features_meta=data.table::fread("2_pops_files/data_from_dropbox/gene_features_metadata.txt")

#Get expression only features
select_cols=colnames(pops_features)[colnames(pops_features) %in% c('ENSGID', pops_features_meta[pops_features_meta$V4=="Expression",]$V1)]
pops_features_expression=pops_features[,..select_cols]
dim(pops_features_expression)
pops_features_expression[1:3, 1:3]
dir.create("2_pops_files/data_from_dropbox/expression_pops/feature_in/", recursive=T)
dir.create("2_pops_files/data_from_dropbox/expression_pops/expression_pops_features_munged/", recursive=T)
data.table::fwrite(pops_features_expression, "2_pops_files/data_from_dropbox/expression_pops/feature_in/expression_only.PoPS.features.txt.gz",
 sep="\t")

#Get Pathway only features
select_cols=colnames(pops_features)[colnames(pops_features) %in% c('ENSGID', pops_features_meta[pops_features_meta$V4=="Pathway",]$V1)]
pops_features_pathway=pops_features[,..select_cols]
dim(pops_features_pathway)
pops_features_pathway[1:3, 1:3]
dir.create("2_pops_files/data_from_dropbox/pathway_pops/feature_in/", recursive=T)
dir.create("2_pops_files/data_from_dropbox/pathway_pops/pathway_pops_features_munged/", recursive=T)
data.table::fwrite(pops_features_pathway, "2_pops_files/data_from_dropbox/pathway_pops/feature_in/pathway_only.PoPS.features.txt.gz", sep="\t")

#Get PPI only features
select_cols=colnames(pops_features)[colnames(pops_features) %in% c('ENSGID', pops_features_meta[pops_features_meta$V4=="PPI",]$V1)]
pops_features_ppi=pops_features[,..select_cols]
dim(pops_features_ppi)
pops_features_ppi[1:3, 1:3]
dir.create("2_pops_files/data_from_dropbox/ppi_pops/feature_in/", recursive=T)
dir.create("2_pops_files/data_from_dropbox/ppi_pops/ppi_pops_features_munged/", recursive=T)
data.table::fwrite(pops_features_ppi, "2_pops_files/data_from_dropbox/ppi_pops/feature_in/ppi_only.PoPS.features.txt.gz", sep="\t")

system("grep -v ppi.control 2_pops_files/data_from_dropbox/full_pops/control.features |grep -v pathways.control > 2_pops_files/data_from_dropbox/expression_pops/control.features")
system("grep ppi.control 2_pops_files/data_from_dropbox/full_pops/control.features > 2_pops_files/data_from_dropbox/ppi_pops/control.features")
system("grep pathways.control 2_pops_files/data_from_dropbox/full_pops/control.features > 2_pops_files/data_from_dropbox/pathway_pops/control.features")
```

#in shell do below

```{bash}
cd /scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir
conda activate GMIP

python ~/GMIP/bin/pops/munge_feature_directory.py \
 --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir 2_pops_files/data_from_dropbox/full_pops/feature_in/ \
 --save_prefix 2_pops_files/data_from_dropbox/full_pops/full_pops_features_munged/pops_features \
 --max_cols 100

python ~/GMIP/bin/pops/munge_feature_directory.py \
 --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir 2_pops_files/data_from_dropbox/expression_pops/feature_in/ \
 --save_prefix 2_pops_files/data_from_dropbox/expression_pops/expression_pops_features_munged/pops_features \
 --max_cols 100

python ~/GMIP/bin/pops/munge_feature_directory.py \
 --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir 2_pops_files/data_from_dropbox/ppi_pops/feature_in/ \
 --save_prefix 2_pops_files/data_from_dropbox/ppi_pops/ppi_pops_features_munged/pops_features \
 --max_cols 100

python ~/GMIP/bin/pops/munge_feature_directory.py \
 --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir 2_pops_files/data_from_dropbox/pathway_pops/feature_in/ \
 --save_prefix 2_pops_files/data_from_dropbox/pathway_pops/pathway_pops_features_munged/pops_features \
 --max_cols 100
```

## PCnet matrix from NAGA framework

The piece of code below is to:\
1. remove all colsums==0 columns from naga matrix.\
2. standardize the columns.\
3. for rows where rowsums==0 add 1 to a control column otherwise 0 (This is how pops does it).\
4. Write the matrix to file.\
5. Write name of column with control feature to control.features file.\

```{r}
library(data.table)
setwd("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir")
features=data.table::fread("2_pops_files/naga/PCnet_matrix_full_for_NAGA.txt.gz")
features[1:10, 1:10]
dim(features)
features_mat=as.matrix(features, rownames=1)
dim(features_mat)
features_mat[1:10, 1:10]
features_mat[is.na(features_mat)]
#Remove columns which are all zeroes
non_zero_columns=colSums(features_mat) != 0
features_mat=features_mat[, non_zero_columns]
dim(features_mat)
features_mat_scaled=scale(features_mat)
features_mat_scaled[is.na(features_mat_scaled)]
features_mat_scaled[1:10, 1:10]

# Add a 'control' column where zero rows have 0 and others have 1
# Identify rows where all columns to be standardized are zero
zero_rows=which(rowSums(features_mat_scaled) == 0)
naga_control_column=rep(1, nrow(features_mat_scaled))
naga_control_column[zero_rows]=0
features_mat_scaled=cbind(naga_control_column, features_mat_scaled)
features_mat_scaled[1:5, 1:5]
features_dt=data.table(ENSGID=rownames(features_mat_scaled), features_mat_scaled)
# Calculate mean and standard deviation for each column
column_means=sapply(features_dt[, -1], mean)
column_sds=sapply(features_dt[, -1], sd)
dim(features_dt)
fwrite(features_dt, "2_pops_files/naga/feature_in/PCnet_from_naga_standardized.tsv.gz", sep = "\t")
sink("2_pops_files/naga/control.features")
cat("naga_control_column")
sink()

```
## GIANT Global top matrix from NetWAS framework

The piece of code below is to:\
1. remove all colsums==0 columns from matrix.\
2. standardize the columns.\
3. for rows where rowsums==0 add 1 to a control column otherwise 0 (This is how pops does it).\
4. Write the matrix to file.\
5. Write name of column with control feature to control.features file.\

```{r}
library(data.table)
setwd("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir")
features=data.table::fread("2_pops_files/netwas/netwas_global_top_matrix_full_for_NAGA.txt.gz")
features[1:10, 1:10]
dim(features)
features_mat=as.matrix(features, rownames=1)
dim(features_mat)
features_mat[1:10, 1:10]
features_mat[is.na(features_mat)]
#Remove columns which are all zeroes
non_zero_columns=colSums(features_mat) != 0
features_mat=features_mat[, non_zero_columns]
dim(features_mat)
features_mat_scaled=scale(features_mat)
features_mat_scaled[is.na(features_mat_scaled)]
features_mat_scaled[1:10, 1:10]

# Add a 'control' column where zero rows have 0 and others have 1
# Identify rows where all columns to be standardized are zero
zero_rows=which(rowSums(features_mat_scaled) == 0)
netwas_control_column=rep(1, nrow(features_mat_scaled))
netwas_control_column[zero_rows]=0
features_mat_scaled=cbind(netwas_control_column, features_mat_scaled)
features_mat_scaled[1:5, 1:5]
features_dt=data.table(ENSGID=rownames(features_mat_scaled), features_mat_scaled)
# Calculate mean and standard deviation for each column
column_means=sapply(features_dt[, -1], mean)
column_sds=sapply(features_dt[, -1], sd)
dim(features_dt)
dir.create("2_pops_files/netwas/feature_in/", showWarnings=T, recursive=T)
fwrite(features_dt, "2_pops_files/netwas/feature_in/GIANT_from_netwas_global_top_standardized.tsv.gz", sep = "\t")
sink("2_pops_files/netwas/control.features")
cat("netwas_control_column")
sink()
```
Preparing for input with pops:
```{bash}
cd /scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir
conda activate GMIP

mkdir -p 2_pops_files/naga/naga_features_munged/
python ~/GMIP/bin/pops/munge_feature_directory.py \
 --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir 2_pops_files/naga/feature_in/ \
 --save_prefix 2_pops_files/naga/naga_features_munged/pops_features \
 --max_cols 100

mkdir -p 2_pops_files/netwas/netwas_features_munged/
python ~/GMIP/bin/pops/munge_feature_directory.py \
 --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir 2_pops_files/netwas/feature_in/ \
 --save_prefix 2_pops_files/netwas/netwas_features_munged/pops_features \
 --max_cols 100
```
