//nextflow run ~/GMIP2/workflows/2_preprocess_wf.nf -with-conda
params.script_dir    = projectDir
params.outdir        = "outdir"
params.prefix        = "out"
params.feature_mat   = null
params.magma_out     = null
params.gene_loc_file = null
params.pval_column   = "P"
params.threshold     = "5e-8"
params.chrom_column  = "CHR"
params.x_id_col      = "ENSGID"
params.y_id_col      = "GENE"
params.folds         = "3"

include { PREPROCESS }  from '../modules/preprocess.nf'

workflow {
  PREPROCESS_otch=PREPROCESS(
    params.script_dir,
    params.pval_column,
    params.threshold,
    params.chrom_column,
    params.x_id_col,
    params.y_id_col,
    params.prefix,
    params.folds,
    params.feature_mat,
    params.magma_out,
    params.gene_loc_file
  )
}
