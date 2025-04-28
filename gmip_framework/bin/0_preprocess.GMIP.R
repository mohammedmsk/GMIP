# Load necessary libraries
library(argparse)
# Parse command-line arguments
parser=ArgumentParser(description="Process X and Y matrices for ML")
parser$add_argument("x_file", help="Path to the X matrix file")
parser$add_argument("y_file", help="Path to the Y matrix file")
parser$add_argument("gene_loc_file", help="Path to the gene_loc file. Made by pipeline in MAGMA directory.")
parser$add_argument("--pval_column", default="P", help="Column name of p-values in Y")
parser$add_argument("--threshold", type="double", default=5e-8, help="P-value threshold for binarization")
parser$add_argument("--chrom_column", default="CHR", help="Column name for chromosomes in Y (used in LOCO)")
parser$add_argument("--x_id_col", default="ENSGID", help="Column name for row IDs in X")
parser$add_argument("--y_id_col", default="GENE", help="Column name for row IDs in Y")
parser$add_argument("--prefix", default="out", help="Prefix for output filenames")
parser$add_argument("--folds", type="integer", default=3, help="Number of folds for Kfold stratification")
parser$add_argument("--nthreads", type="integer", default=4, help="Number of threads")
args=parser$parse_args()

library(data.table)
library(caret)
library(parallel)

# Function to subset X and Y to have the same rows and adjust order
subset_XY=function(x, y, x_id_col="ENSGID", y_id_col="GENE") {
  common_ids=intersect(x[[x_id_col]], y[[y_id_col]])
  if (length(common_ids) < 100) {
    stop("Not enough common genes between X and Y, did you specify ID column names correctly?")
  }
  # Subset and reorder X and Y based on common_ids
  Xtest2=x[!get(x_id_col) %in% common_ids]
  x=x[get(x_id_col) %in% common_ids]
  y=y[get(y_id_col) %in% common_ids]
  # Reorder rows to match the order of common_ids
  x=x[match(common_ids, x[[x_id_col]])]
  y=y[match(common_ids, y[[y_id_col]])]
  if (!identical(x$ENSGID, y$GENE)) {
    stop("X and Y IDs do not match.")
  }
  print(paste0("Dimensions of X: ", dim(x)))
  print(paste0("Dimensions of Y: ", dim(y)))
  print(paste0("Dimensions of Xtest2: ", dim(Xtest2)))
  return(list(x=x, y=y, Xtest2=Xtest2))
}

# Function to perform stratified K-fold split
perform_stratified_split=function(x, y, Xtest2, folds=3, prefix="out", nthreads=4) {
  set.seed(786786)
  fold_indices=createFolds(y$binarized_P, k=folds, list=TRUE, returnTrain=TRUE)
  output_data_list=mclapply(seq_along(fold_indices), function(i) {
    Xtrain=x[fold_indices[[i]], ]
    Ytrain=y[fold_indices[[i]], ]
    Xtest=x[-fold_indices[[i]], ]
    Ytest=y[-fold_indices[[i]], ]
    fold_suffix=paste0("_fold", i)
    Xtrain_filename=paste0(prefix, "_Xtrain", fold_suffix, ".tsv")
    Ytrain_filename=paste0(prefix, "_Ytrain", fold_suffix, ".tsv")
    Xtest_filename=paste0(prefix, "_Xtest", fold_suffix, ".tsv")
    Ytest_filename=paste0(prefix, "_Ytest", fold_suffix, ".tsv")
    Xtest2_filename=paste0(prefix, "_Xtest2", fold_suffix, ".tsv")
    fwrite(Xtrain, Xtrain_filename, sep="\t", nThread=2)
    fwrite(Ytrain, Ytrain_filename, sep="\t", nThread=2)
    fwrite(Xtest, Xtest_filename, sep="\t", nThread=2)
    fwrite(Ytest, Ytest_filename, sep="\t", nThread=2)
    fwrite(Xtest2, Xtest2_filename, sep="\t", nThread=2)
    data.table(fold=paste0("fold", i), XtrainFileName=Xtrain_filename, XtestFileName=Xtest_filename,
      YtrainFileName=Ytrain_filename, YtestFileName=Ytest_filename, Xtest2FileName=Xtest2_filename)
  }, mc.cores=nthreads)

  fwrite(rbindlist(output_data_list, use.names=T, fill=T), paste0(prefix, "_foldFiles.csv"), sep=",")
  return("Done!")
}

# Function to perform leave-one-chromosome-out (LOCO) split
perform_loco_split=function(x, y, Xtest2, gene_loc, chrom_column="CHR", prefix="out", nthreads=4) {
  unique_chroms=unique(y[[chrom_column]])
  fwrite(data.table(Chroms=unique_chroms), paste0(prefix, "_chroms.txt"), col.names=F)

  output_data_list=mclapply(unique_chroms, function(chrom) {
    Xtrain=x[y[[chrom_column]] != chrom, ]
    Ytrain=y[y[[chrom_column]] != chrom, ]
    Xtest=x[y[[chrom_column]] == chrom, ]
    Ytest=y[y[[chrom_column]] == chrom, ]
    og_cols_Xtest2=colnames(Xtest2)
    Xtest2=merge(gene_loc, Xtest2, by=colnames(Xtest2)[1])
    Xtest2=Xtest2[Xtest2[[chrom_column]] == chrom, ]
    Xtest2=Xtest2[, ..og_cols_Xtest2]
    chrom_suffix=paste0("_chrom", chrom)
    Xtrain_filename=paste0(prefix, "_Xtrain", chrom_suffix, ".tsv")
    Ytrain_filename=paste0(prefix, "_Ytrain", chrom_suffix, ".tsv")
    Xtest_filename=paste0(prefix, "_Xtest", chrom_suffix, ".tsv")
    Ytest_filename=paste0(prefix, "_Ytest", chrom_suffix, ".tsv")
    Xtest2_filename=paste0(prefix, "_Xtest2", chrom_suffix, ".tsv")
    fwrite(Xtrain, Xtrain_filename, sep="\t")
    fwrite(Ytrain, Ytrain_filename, sep="\t")
    fwrite(Xtest, Xtest_filename, sep="\t")
    fwrite(Ytest, Ytest_filename, sep="\t")
    fwrite(Xtest2, Xtest2_filename, sep="\t")
    data.table(fold=paste0("chrom", chrom), XtrainFileName=Xtrain_filename, XtestFileName=Xtest_filename,
      YtrainFileName=Ytrain_filename, YtestFileName=Ytest_filename, Xtest2FileName=Xtest2_filename)
  }, mc.cores=nthreads)

  fwrite(rbindlist(output_data_list, use.names=T, fill=T), paste0(prefix, "_chromFiles.csv"), sep=",")
  return("Done!")
}

# Main function
process_data=function(x_file, y_file, gene_loc_file, pval_column="P", threshold=5e-8, chrom_column="CHR", x_id_col="ENSGID",
  y_id_col="GENE", prefix="out", folds=3, nthreads=4) {

  # Read X and Y matrices
  x=fread(x_file)
  y=fread(y_file)
  gene_loc=fread(gene_loc_file)

  # Subset X and Y to have the same rows
  subsetted_data=subset_XY(x, y, x_id_col, y_id_col)
  x=subsetted_data$x
  y=subsetted_data$y
  Xtest2=subsetted_data$Xtest2

  print(x[1:5, 1:5])

  # Process Y matrix to add neg_log_pval and binary columns
  y[, neg_log_P := -log10(get(pval_column))]
  y[, binarized_P := ifelse(neg_log_P > -log10(threshold), 1, 0)]
  print(y[1:5,])

  # These steps are to make pops work smoothly
  y$ENSGID=y[[y_id_col]]
  if ("ZSTAT" %in% colnames(y)) {
    y$Score=y$ZSTAT
  } else {
    y$Score=y$neg_log_P
  }

  #X and Y will be same for test and train when noCV
  Xtrain_filename=paste0(prefix, "_Xtrain_full.tsv")
  Ytrain_filename=paste0(prefix, "_Ytrain_full.tsv")
  Xtest2_filename=paste0(prefix, "_Xtest2_full.tsv")
  fwrite(x, Xtrain_filename, sep="\t")
  fwrite(y, Ytrain_filename, sep="\t")
  fwrite(Xtest2, Xtest2_filename, sep="\t")
  dt=data.table(fold="full", XtrainFileName=Xtrain_filename, XtestFileName=Xtrain_filename,
    YtrainFileName=Ytrain_filename, YtestFileName=Ytrain_filename, Xtest2FileName=Xtest2_filename)
  fwrite(dt, paste0(prefix, "_noCV_Files.csv"), sep=",")

  # Perform stratified K-fold split and write files
  perform_stratified_split(x, y, Xtest2, folds=folds, prefix=prefix, nthreads=nthreads)

  # Perform leave-one-chromosome-out (LOCO) split and write files
  perform_loco_split(x, y, Xtest2, gene_loc, chrom_column=chrom_column, prefix=prefix, nthreads=nthreads)

  return("Done making preprocessing files.")
}

# args=list(x_file="/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/data_from_dropbox/ppi_pops/feature_in/ppi_only.PoPS.features.txt",
#   y_file="/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/SCZ.magma.genes.out", gene_loc_file="/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/gene_loc.txt", pval_column="P", threshold=5e-8, chrom_column="CHR", x_id_col="ENSGID", y_id_col="GENE", prefix="out", folds=3)

# Run the main function with parsed arguments
process_data(args$x_file, args$y_file, args$gene_loc_file, args$pval_column, args$threshold, args$chrom_column,
  args$x_id_col, args$y_id_col, args$prefix, args$folds, args$nthreads)


