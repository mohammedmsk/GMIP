process EVAL {
  publishDir "${params.outdir}/6_eval_results/${outprefix}/", mode: 'copy', overwrite: true

  input:
    path(script_dir)
    path(ytestFile)
    path(ypredFile)
    val(outprefix)
    path(gold_standardFile)
    val(GENE_COL1)
    val(SCORE_COL1)
    val(BINARY_COL1)
    val(GENE_COL2)
    val(SCORE_COL2)

  output:
    path("${outprefix}*")

  script:
  """
  Rscript ${script_dir}/bin/eval.R \\
    --gene_col1 ${GENE_COL1} \\
    --score_col1 ${SCORE_COL1} \\
    --binary_col1 ${BINARY_COL1} \\
    --gene_col2 ${GENE_COL2} \\
    --score_col2 ${SCORE_COL2} \\
    ${ytestFile} \\
    ${ypredFile} \\
    ${outprefix} \\
    ${gold_standardFile}
  """
}
