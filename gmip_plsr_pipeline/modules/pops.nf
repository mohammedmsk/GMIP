process POPS {
  publishDir "${params.pops.outdir}", mode: 'copy', overwrite: true

  input:
    tuple(
      path(script_dir),
      path(gene_annot_path),
      path(feature_mat_pathdir),
      val(feature_mat_prefix),
      path(control_features_path),
      val(chr)
    )
    path(magma_out)
    path(magma_raw)

  output:
    path("${params.pops.prefix}.test_chr${chr}.coefs",     emit: pops_coefs)
    path("${params.pops.prefix}.test_chr${chr}.genes",     emit: pops_genes)
    path("${params.pops.prefix}.test_chr${chr}.marginals", emit: pops_marginals)
    path("${params.pops.prefix}.test_chr${chr}.matdata",   emit: pops_matdata)
    path("${params.pops.prefix}.test_chr${chr}.preds",     emit: pops_preds)
    path("${params.pops.prefix}.test_chr${chr}.preds2",    emit: pops_preds2)
    path("${params.pops.prefix}.test_chr${chr}.traindata", emit: pops_traindata)
    path("${params.pops.prefix}.test_chr${chr}.log",       emit: pops_log)

  script:
  def magma_prefix = magma_out.toString().replaceAll(/\.genes\.out$/, "")
  """
  filteredChromosomes=\$(echo "${params.chromosomes}" | tr ' ' '\n' | grep -wv "^${chr}" | tr '\n' ' ')

  echo "Testing for chromosome: chr${chr}"

  pops_options=''

  if [ "${params.pops.use_magma_covariates}" = "true" ] && [ "${params.collapse_loci}" = "false" ]; then
    pops_options="\${pops_options} --use_magma_covariates"
  else
    pops_options="\${pops_options} --ignore_magma_covariates"
  fi

  if [ "${params.pops.use_magma_error_cov}" = "true" ] && [ "${params.collapse_loci}" = "false" ]; then
    pops_options="\${pops_options} --use_magma_error_cov"
  else
    pops_options="\${pops_options} --ignore_magma_error_cov"
  fi

  if [ "${params.pops.project_out_covariates_remove_hla}" = "true" ]; then
    pops_options="\${pops_options} --project_out_covariates_remove_hla"
  else
    pops_options="\${pops_options} --project_out_covariates_keep_hla"
  fi

  if [ "${params.pops.feature_selection_remove_hla}" = "true" ]; then
    pops_options="\${pops_options} --feature_selection_remove_hla"
  else
    pops_options="\${pops_options} --feature_selection_keep_hla"
  fi

  if [ "${params.pops.training_remove_hla}" = "true" ]; then
    pops_options="\${pops_options} --training_remove_hla"
  else
    pops_options="\${pops_options} --training_keep_hla"
  fi

  if [ ! -z "${control_features_path}" ]; then
    pops_options="\${pops_options} --control_features_path ${control_features_path}"
  else
    pops_options="\${pops_options}"
  fi

  num_feature_chunks=`ls ${feature_mat_pathdir}/${feature_mat_prefix}.mat.*.npy|wc -l`

  which python

  python ${script_dir}/bin/pops/pops.py \
    --gene_annot_path ${gene_annot_path} \
    --feature_mat_prefix ${feature_mat_pathdir}/${feature_mat_prefix} \
    --num_feature_chunks \${num_feature_chunks} \
    --magma_prefix ${magma_prefix} \
    --feature_selection_p_cutoff 0.05 \
    --method ${params.pops.method} \
    --out_prefix ${params.pops.prefix}.test_chr${chr} \
    --verbose --random_seed 42 --save_matrix_files \
    --training_chromosomes \${filteredChromosomes} \
    --project_out_covariates_chromosomes \${filteredChromosomes} \
    --feature_selection_chromosomes \${filteredChromosomes} \
    \${pops_options} &> ${params.pops.prefix}.test_chr${chr}.log


  awk -v chr=${chr} '{if(\$3==chr) print \$1}' ${gene_annot_path} > ${params.pops.prefix}.test_chr${chr}.genes
  grep -w -Ff ${params.pops.prefix}.test_chr${chr}.genes ${params.pops.prefix}.test_chr${chr}.preds|\
    sort -k2,2gr > ${params.pops.prefix}.test_chr${chr}.preds2

  """
}
