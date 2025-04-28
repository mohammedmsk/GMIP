process GMIP {
  publishDir "${params.gmip.outdir}", mode: 'copy', overwrite: true

  input:
    path(script_dir)
    path(gene_annot_path)
    path(trainFile)
    path(matFile)
    path(marginalFile)
    val(pprefix)
    val(method)
    val(noFeatureClustering)

  output:
    path("${pprefix}*.test_chr*.coefs",  emit: gmip_coefs)
    path("${pprefix}*.test_chr*.genes",  emit: gmip_genes)
    path("${pprefix}*.test_chr*.preds",  emit: gmip_preds)
    path("${pprefix}*.test_chr*.preds2", emit: gmip_preds2)
    path("${pprefix}*.test_chr*.log",    emit: gmip_log)

  script:
  def chr=trainFile.toString().replaceAll(/^.*\.test_chr(\d+)\.traindata$/, '$1')
  """
  if [ "${noFeatureClustering}" = "true" ]; then
    gmip_options=" --noFeatureClustering"
  else
    gmip_options=''
  fi

  python ${script_dir}/bin/gmip.py \
    --trainFile ${trainFile} \
    --matFile ${matFile} \
    --marginalFile ${marginalFile} \
    --outPrefix ${pprefix}.test_chr${chr} \
    --methodName ${method} \
    \${gmip_options} &> ${pprefix}.test_chr${chr}.log

  awk -v chr=${chr} '{if(\$3==chr) print \$1}' ${gene_annot_path} > ${pprefix}.test_chr${chr}.genes
  grep -w -Ff ${pprefix}.test_chr${chr}.genes ${pprefix}.test_chr${chr}.preds|\
    sort -k2,2gr > ${pprefix}.test_chr${chr}.preds2
  """
}

