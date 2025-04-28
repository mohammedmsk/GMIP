process COLLAPSE_GENES_TO_LOCI {
  publishDir "${params.magma.outdir}", mode: 'copy', overwrite: true

  input:
      path(yfile)
      val(id)
      val(chr)
      val(start)
      val(end)
      val(score)
      val(wd_sz)

  output:
    path("${yfile}.collapsed.genes.out",  emit: magma_out)
    path("${yfile}.collapsed.genes.raw",  emit: magma_raw)
    path("${yfile}.loci",                 emit: magma_loci)

  script:
  """
  python ${params.script_dir}/bin/collapse_genes_to_loci.py \\
    --yfile ${yfile} \\
    --id ${id} \\
    --chr ${chr} \\
    --start ${start} \\
    --end ${end} \\
    --score ${score} \\
    --wd_sz ${wd_sz}

  #Making an empty magma_raw file as it is needed for pops, not actually used at all
  cat > ${yfile}.collapsed.genes.raw
  """
}
