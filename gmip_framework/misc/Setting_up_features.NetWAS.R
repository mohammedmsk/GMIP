make_full=function(x){
  my_cols=colnames(x)
  my_rows=rownames(x)
  ##Add missing rows to the matrix
  to_add_as_rows=setdiff(my_cols, my_rows)
  if(length(to_add_as_rows)!=0){
    m1=c(rep.int(0, length(to_add_as_rows)*length(my_cols)))
    dim(m1)=c(length(to_add_as_rows), length(my_cols))
    colnames(m1)=colnames(x)
    rownames(m1)=to_add_as_rows
    x=rbind(x, m1)
  }
  ##Add missing cols to the matrix
  to_add_as_cols=setdiff(my_rows, my_cols)
  if(length(to_add_as_cols)!=0){
    m2=rep.int(0, nrow(x)*length(to_add_as_cols))
    dim(m2)=c(nrow(x), length(to_add_as_cols))
    colnames(m2)=to_add_as_cols
    x=cbind(x, m2)
  }
  return(x)
}

library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(VennDiagram)
library(gridExtra)
library(ndexr)
library(igraph)
library(RCy3)

setwd("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files")
dir.create(path="netwas", showWarnings=F, recursive=T)
setwd("netwas")
genes_from_pops=fread("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir/static_gits/pops/example/data/utils/gene_annot_jun10.txt")
genes_to_keep=genes_from_pops$ENSGID

#This fro HumanBase download. Global network for netwas. Top edges
file_url="https://s3-us-west-2.amazonaws.com/humanbase/networks/global_top.gz"

# Download the file
download.file(file_url, "global_top.gz", method="auto")

#Read in the files using data.tables
elist=data.table::fread("global_top.gz", header=F)
colnames(elist)=c("Target", "Source", "Weight")
elist$Target=as.character(elist$Target)
elist$Source=as.character(elist$Source)
#This command converts the elist into matrix-like data. Fastest in R so stick with it.
x_adj=tidyr::spread(elist, key=Source, value=Weight, fill=0, convert=T)
x_adj[1:5,1:5]
dim(x_adj)
colnames(x_adj)[1]="ID"
#Now matching row
data.table::setorder(x_adj, ID)
data.table::setcolorder(x_adj, c("ID",base::sort(colnames(x_adj)[-1])))
x_adj[1:10, 1:5]
dim(x_adj)
x_adj_mat=as.matrix(x_adj[,-1])
rownames(x_adj_mat)=x_adj$ID
x_adj_full=make_full(x_adj_mat)
dim(x_adj_full)
x_adj_full[1:10, 1:10]
x_adj_full=x_adj_full[,rownames(x_adj_full)]
identical(rownames(x_adj_full), colnames(x_adj_full))

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve gene symbols
gene_info <- getBM(
  filters = "entrezgene_id",
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  values = rownames(x_adj_full),
  mart = ensembl
)

dim(gene_info)
length(unique(rownames(x_adj_full)))
x_adj_m=merge(gene_info, x_adj_full, by.x="entrezgene_id", by.y=0)
dim(x_adj_m)
x_adj_m[1:10, 1:10]
fwrite(x_adj_m, "global_top_matrix.txt.gz", sep="\t")
x_adj_m2=x_adj_m[,-1]
x_adj_m3=as.matrix(x_adj_m2[,-1])
rownames(x_adj_m3)=x_adj_m2[,1]
x_adj_m3[1:5, 1:5]
x_adj_m3=t(x_adj_m3)
x_adj_m3[1:5, 1:5]

# Retrieve gene symbols
gene_info2 <- getBM(
  filters = "entrezgene_id",
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  values = rownames(x_adj_m3),
  mart = ensembl
)

dim(gene_info2)
identical(gene_info, gene_info2)

length(unique(rownames(x_adj_m3)))
x_adj_mm=merge(gene_info2, x_adj_m3, by.x="entrezgene_id", by.y=0)
dim(x_adj_mm)
x_adj_mm[1:10, 1:10]
x_adj_mm2=x_adj_mm[,-1]
x_adj_mm3=as.matrix(x_adj_mm2[,-1])
rownames(x_adj_mm3)=x_adj_mm2[,1]
x_adj_mm3[1:5, 1:5]
dim(x_adj_mm3)
identical(rownames(x_adj_mm3), colnames(x_adj_mm3))
x_adj_mm3=x_adj_mm3[,rownames(x_adj_mm3)]
identical(rownames(x_adj_mm3), colnames(x_adj_mm3))

#ndex_conn <- ndex_connect()
#network_id <- "f93f402c-86d4-11e7-a10d-0ac135e8bacf"  # Replace with your network ID
#cx_network <- ndex_get_network(ndex_conn, network_id)

venn.plot <- venn.diagram(
  x = list(
    pops_genes = genes_to_keep,
    netwas_genes = rownames(x_adj_mm3)
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "blue")
)

grid.draw(venn.plot)

final_genes=intersect(genes_to_keep, rownames(x_adj_mm3))


# Identify rows with non-zero sum
non_zero_rows <- rowSums(x_adj_mm3) != 0

# Identify columns with non-zero sum
non_zero_cols <- colSums(x_adj_mm3) != 0

fwrite(data.table(ID=rownames(x_adj_mm3), x_adj_mm3), "netwas_global_top_matrix_full_for_NAGA.txt.gz", sep="\t")

# Subset the matrix
filtered_mat <- x_adj_mm3[non_zero_rows, non_zero_cols]
dim(x_adj_mm3)
dim(filtered_mat)

filtered_mat_dt=data.table(ID=rownames(filtered_mat))
filtered_mat_dt=cbind(filtered_mat_dt, as.data.table(filtered_mat))
filtered_mat_dt[1:5, 1:5]
dim(filtered_mat_dt)
cols_to_standardize=names(filtered_mat_dt)[-1]

# Step 4: Standardize the columns except the first one using scale()
filtered_mat_dt[, (cols_to_standardize) := lapply(.SD, scale), .SDcols=cols_to_standardize]


