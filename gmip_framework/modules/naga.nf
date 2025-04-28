process NAGA {
  publishDir "${params.outdir}/3_naga/${prefix}", mode: 'copy', overwrite: true

  input:
    path(script_dir)
    path(feature_network)
    path(preprocess_files)
    val(prefix)
    tuple(val(fold), val(XtrainFileName), val(XtestFileName), val(YtrainFileName), val(YtestFileName), val(Xtest2FileName))
    val(ScoreColname)
    val(IDColname)

  output:
    path("${prefix}.*.genes.nagaRes.tsv", emit: naga_res)
    path("${prefix}.*.test_genes.list",   emit: naga_test_genes_list)
    path("${prefix}.*.yhat",              emit: naga_yhat)

  script:
  """
    python ${script_dir}/bin/run_naga.py \
      --YtrainFileName ${YtrainFileName} \
      --XnetworkFileName ${feature_network} \
      --ScoreColname ${ScoreColname} \
      --IDColname ${IDColname} \
      --Prefix ${prefix}.\$(basename ${YtrainFileName} .tsv)
    
    echo ENSGID > ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list
    cat ${YtestFileName}|cut -f1|tail -n+2 >> ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list
    cat ${Xtest2FileName}|cut -f1|tail -n+2 >> ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list

    grep -w -Ff ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list\
      ${prefix}.\$(basename ${YtrainFileName} .tsv).genes.nagaRes.lim.tsv\
      > ${prefix}.\$(basename ${YtrainFileName} .tsv).yhat
  """
}
