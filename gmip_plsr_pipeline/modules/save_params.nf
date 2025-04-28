process saveParams {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  output:
  path("${params.prefix}.params.txt")

  script:
  def paramsList = params.collect { k, v -> "$k : $v" }
  def paramsString = paramsList.join('\n')
  """
  echo "Workflow Parameters" > ${params.prefix}.params.txt
  echo "-------------------" >> ${params.prefix}.params.txt
  echo "" >> ${params.prefix}.params.txt
  echo "${paramsString}" >> ${params.prefix}.params.txt
  """
}
