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

setwd("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files")
dir.create(path="naga", showWarnings=F, recursive=T)
setwd("naga")
genes_from_pops=as.data.frame(fread("/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir/static_gits/pops/example/data/utils/gene_annot_jun10.extended.txt"))

#Read in the files using data.tables
x_adj=data.table::fread("csr_matrix_with_default_labels.txt.gz")
x_adj[1:5,1:5]
dim(x_adj)
colnames(x_adj)[1]="ID"
#Now matching row
data.table::setorder(x_adj, ID)
data.table::setcolorder(x_adj, c("ID",base::sort(colnames(x_adj)[-1])))
x_adj[1:10, 1:5]
dim(x_adj)
#convert to matrix
x_adj_mat=as.matrix(x_adj[,-1])
rownames(x_adj_mat)=x_adj$ID
dim(x_adj_mat)
x_adj_full=make_full(x_adj_mat)
dim(x_adj_full)
x_adj_full[1:10, 1:10]
head(genes_from_pops)
x_adj2=merge(genes_from_pops, x_adj_full, by.x="NAME", by.y=0)
x_adj2[1:10, 1:10]
dim(x_adj2)
x_adj3=x_adj2[,c(colnames(genes_from_pops), x_adj2$NAME)]
dim(x_adj3)
rownames(x_adj3)=x_adj3$ENSGID
cols_to_rv=which(colnames(genes_from_pops) %in% colnames(x_adj3))
x_adj3=x_adj3[,-cols_to_rv]
x_adj3[1:10, 1:10]
dim(x_adj3)
colnames(x_adj3)=rownames(x_adj3)
x_adj4=merge(genes_from_pops, x_adj3, by.x="ENSGID", by.y=0, all.x=T)
x_adj4[1:5, 1:15]
x_adj4[is.na(x_adj4)]=0
rownames(x_adj4)=x_adj4$ENSGID
cols_to_rv=which(colnames(genes_from_pops) %in% colnames(x_adj4))
x_adj4=x_adj4[,-cols_to_rv]
x_adj4[1:10, 1:10]
dim(x_adj4)
x_adj4_mat=as.matrix(x_adj4)
x_adj4_full=make_full(x_adj4_mat)
dim(x_adj4_full)
x_adj4_full[1:10, 1:10]
x_adj4_full=x_adj4_full[,rownames(x_adj4_full)]
identical(rownames(x_adj4_full), colnames(x_adj4_full))
x_adj4_full[1:10, 1:10]
fwrite(data.table(ENSGID=rownames(x_adj4_full), x_adj4_full), "PCnet_matrix_full_for_NAGA.txt.gz", sep="\t")
