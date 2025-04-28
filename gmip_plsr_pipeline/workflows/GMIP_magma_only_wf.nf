//
include { MAGMA } from '../modules/magma.nf'
workflow {

  MAGMA_otch=MAGMA(
    params.magma.executable,
    params.magma_input,
    params.magma.ref_pathdir,
    params.magma.bfile_prefix,
    params.gene_annot_path,
    params.magma.outdir,
    params.magma.prefix
  )
}
