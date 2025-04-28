//nextflow run ~/GMIP2/workflows/3_pops_wf.nf -with-conda
params.base_refdir             = "/scratch/04179/mshabb/04488464786465dwkljhefkhbefjuh/refdir"
params.script_dir              = "/home1/04179/mshabb/GMIP2"
params.conda_path              = "/work2/04179/mshabb/ls6/miniforge3/envs/"
params.gene_annot_path         = "${params.script_dir}/bin/pops/example/data/utils/gene_annot_jun10.txt"
params.feature_mat_pathdir     = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/ppi_pops_features_munged"
params.feature_mat_prefix      = "pops_features"
params.control_features_path   = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/control.features"
params.outdir                  = "outdir"
params.prefix                  = "NAFLD_popsFeatures.pops"
params.y_train                 = "/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/out_Ytrain_chrom10.tsv.gz"
params.y_test                  = "/scratch/04179/mshabb/GMIP2/data/MAGMA/SCZ/outdir/1_magma/out_Ytest_chrom10.tsv.gz"

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
