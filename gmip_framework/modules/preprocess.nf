process PREPROCESS {
  publishDir "${params.outdir}/2_preprocess", mode: 'copy', overwrite: true

  input:
    path(script_dir)
    val(pval_column)
    val(threshold)
    val(chrom_column)
    val(x_id_col)
    val(y_id_col)
    val(prefix)
    val(folds)
    path(feature_mat)
    path(magma_out)
    path(gene_loc_file)

  output:
    path("${prefix}_chroms.txt",         emit:preprocess_chroms)
    path("${prefix}_chromFiles.csv",     emit:preprocess_chrom_csv)
    path("${prefix}_*_chrom*.tsv*",      emit:preprocess_chrom_files)
    path("${prefix}_foldFiles.csv",      emit:preprocess_fold_csv)
    path("${prefix}_*_fold*.tsv*",       emit:preprocess_fold_files)
    path("${prefix}_noCV_Files.csv",     emit:preprocess_nocv_csv)
    path("${prefix}*_full.tsv*",         emit:preprocess_nocv_files)
    path("${prefix}*_Ytrain_full.tsv*",  emit:preprocess_Ytrain_full)

  script:
  """
  Rscript ${script_dir}/bin/0_preprocess.GMIP.R \\
    --pval_column ${pval_column} \\
    --threshold ${threshold} \\
    --chrom_column ${chrom_column} \\
    --x_id_col ${x_id_col} \\
    --y_id_col ${y_id_col}\\
    --prefix ${prefix} \\
    --folds ${folds} \\
    ${feature_mat} \\
    ${magma_out} \\
    ${gene_loc_file} \\
    --nthreads ${task.cpus}
  """
}
