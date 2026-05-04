//nextflow run ~/GMIP2/workflows/3_pops_wf.nf -with-conda
params.base_refdir             = "${projectDir}/refdir"
params.script_dir              = projectDir
params.gene_annot_path         = "${projectDir}/bin/pops/example/data/utils/gene_annot_jun10.txt"
params.feature_mat_pathdir     = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/ppi_pops_features_munged"
params.feature_mat_prefix      = "pops_features"
params.control_features_path   = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/control.features"
params.outdir                  = "outdir"
params.prefix                  = "NAFLD_popsFeatures.pops"
params.y_train                 = null
params.y_test                  = null

include { POPS }  from '../modules/pops.nf'

workflow {
  POPS_otch=POPS(
    params.script_dir,
    params.gene_annot_path,
    params.feature_mat_pathdir,
    params.feature_mat_prefix,
    params.control_features_path,
    params.y_train, params.y_test
  )
}
