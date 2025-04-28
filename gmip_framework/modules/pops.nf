process POPS {
  publishDir "${params.outdir}/3_pops/${prefix}", mode: 'copy', overwrite: true

  input:
    path(script_dir)
    path(gene_annot_path)
    path(feature_mat_pathdir)
    val(feature_mat_prefix)
    path(control_features_path)
    path(preprocess_files)
    val(prefix)
    tuple(val(fold), val(XtrainFileName), val(XtestFileName), val(YtrainFileName), val(YtestFileName), val(Xtest2FileName))

  output:
    path("${prefix}.*.coefs",           emit: pops_coefs)
    path("${prefix}.*.log",             emit: pops_log)
    path("${prefix}.*.marginals",       emit: pops_marginals)
    path("${prefix}.*.matdata",         emit: pops_matdata)
    path("${prefix}.*.preds",           emit: pops_preds)
    path("${prefix}.*.test_genes.list", emit: pops_test_genes_list)
    path("${prefix}.*.traindata",       emit: pops_traindata)
    path("${prefix}.*.yhat",            emit: pops_yhat)

  script:
  """
    num_feature_chunks=`ls ${feature_mat_pathdir}/${feature_mat_prefix}.mat.*.npy|wc -l`

    python ${script_dir}/bin/pops/pops.py \
      --gene_annot_path ${gene_annot_path} \
      --feature_mat_prefix ${feature_mat_pathdir}/${feature_mat_prefix} \
      --num_feature_chunks \${num_feature_chunks} \
      --feature_selection_p_cutoff 0.05 \
      --method ridge \
      --out_prefix ${prefix}.\$(basename ${YtrainFileName} .tsv) \
      --verbose --random_seed 42 --save_matrix_files \
      --y_path ${YtrainFileName} \
      --control_features_path ${control_features_path} \
      &> ${prefix}.\$(basename ${YtrainFileName} .tsv).log

    echo ENSGID > ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list
    cat ${YtestFileName}|cut -f1|tail -n+2 >> ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list
    cat ${Xtest2FileName}|cut -f1|tail -n+2 >> ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list

    grep -w -Ff ${prefix}.\$(basename ${YtrainFileName} .tsv).test_genes.list\
      ${prefix}.\$(basename ${YtrainFileName} .tsv).preds |sed 's/PoPS_Score/PredictedScore/'\
      > ${prefix}.\$(basename ${YtrainFileName} .tsv).yhat
  """
}
