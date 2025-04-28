//nextflow run ~/GMIP2/workflows/2_3_preprocess_pops_wf.nf -with-conda

params.base_refdir  = "/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir"
params.conda_path   = "/work2/04179/mshabb/ls6/miniforge3/envs/"
params.script_dir   = "/home1/04179/mshabb/GMIP2"
params.outdir       = "outdir"
params.prefix       = "out"

//preprocess related args
params.feature_mat   = "/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/data_from_dropbox/ppi_pops/feature_in/ppi_only.PoPS.features.txt.gz"
params.magma_out     = "/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/SCZ.magma.genes.out"
params.gene_loc_file = "/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/gene_loc.txt"
params.pval_column   = "P"
params.threshold     = "5e-8"
params.chrom_column  = "CHR"
params.x_id_col      = "ENSGID"
params.y_id_col      = "GENE"
params.folds         = "3"

//pops related args
params.gene_annot_path         = "${params.script_dir}/bin/pops/example/data/utils/gene_annot_jun10.txt"
params.feature_mat_pathdir     = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/ppi_pops_features_munged"
params.feature_mat_prefix      = "pops_features"
params.control_features_path   = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/control.features"

include { PREPROCESS }  from '../modules/preprocess.nf'
include { POPS }  from '../modules/pops.nf'

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

  def lineCount=0
  ml_inch=PREPROCESS_otch.preprocess_chrom_csv.splitText().map { line ->
    if (lineCount++ == 0) {
      return null // Skip the header
    } else {
      def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
      return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName)
    }
  }

  POPS(
    params.script_dir,
    params.gene_annot_path,
    params.feature_mat_pathdir,
    params.feature_mat_prefix,
    params.control_features_path,
    PREPROCESS_otch.preprocess_chrom_files.collect(),
    ml_inch
  )
}
