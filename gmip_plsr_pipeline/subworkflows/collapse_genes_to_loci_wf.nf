include { COLLAPSE_GENES_TO_LOCI } from "../modules/collapse_genes_to_loci.nf"

workflow {
  main:
  COLLAPSE_GENES_TO_LOCI(
    params.yfile,
    params.collapse_genes_to_loci.id,
    params.collapse_genes_to_loci.chr,
    params.collapse_genes_to_loci.start,
    params.collapse_genes_to_loci.end,
    params.collapse_genes_to_loci.score,
    params.collapse_genes_to_loci.wd_sz
  )
}
