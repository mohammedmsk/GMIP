//nextflow run ~/GMIP2/workflows/2_preprocess_wf.nf -with-conda
params.conda_path    = "/work2/04179/mshabb/ls6/miniforge3/envs/"
params.script_dir    = "/home1/04179/mshabb/GMIP2"
params.outdir        = "outdir"
params.prefix        = "out"
params.feature_mat   = "/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/data_from_dropbox/ppi_pops/feature_in/ppi_only.PoPS.features.txt.gz"
params.magma_out     = "/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/SCZ.magma.genes.out"
params.gene_loc_file = "/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/gene_loc.txt"
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
