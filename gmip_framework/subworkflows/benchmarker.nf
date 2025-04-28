workflow BENCHMARKER {

  take:
  genes_top
  pprefix
  label

  main:
  // Define a channel with all autosomal chromosome numbers
  chromosome_ch = Channel.from(1..22)

  // Call BENCHMARKER_PART1 for each chromosome
  BENCHMARKER_PART1_inch = chromosome_ch.map { chr ->
    tuple(
      params.script_dir,
      params.gene_annot_path,
      params.benchmarker.baseline_pathdir,
      params.benchmarker.baseline_prefix,
      params.benchmarker.bfile_pathdir,
      params.benchmarker.bfile_prefix,
      params.benchmarker.print_snps,
      params.benchmarker.window_size_in_kb,
      chr)
  }

  BENCHMARKER_PART1_otch = BENCHMARKER_PART1(BENCHMARKER_PART1_inch, genes_top, pprefix, label)

  BENCHMARKER_PART2(
    BENCHMARKER_PART1_otch.for_bm_part2_annot.collect(),
    BENCHMARKER_PART1_otch.for_bm_part2_ld.collect(),
    BENCHMARKER_PART1_otch.for_bm_part2_m.collect(),
    BENCHMARKER_PART1_otch.for_bm_part2_m_5_50.collect(),
    BENCHMARKER_PART1_otch.for_bm_part2_log.collect(),
    params.script_dir,
    params.munged_gwas,
    params.benchmarker.baseline_pathdir,
    params.benchmarker.baseline_prefix,
    params.benchmarker.gene_window_pathdir,
    params.benchmarker.gene_window_prefix,
    params.benchmarker.weights_pathdir,
    params.benchmarker.weights_prefix,
    params.benchmarker.frq_pathdir,
    params.benchmarker.frq_prefix,
    params.benchmarker.window_size_in_kb,
    pprefix,
    label
  )
}

process BENCHMARKER_PART1 {
  publishDir "${params.outdir}/5_bm_results/${pprefix}/part_1/", mode: 'copy', overwrite: true

  input:
    tuple(
      path(script_dir),
      path(gene_annot_path),
      path(baseline_pathdir),
      val(baseline_prefix),
      path(bfile_pathdir),
      val(bfile_prefix),
      path(print_snps),
      val(window_size),
      val(chr)
    )
    path(prio_gene_list)
    val(pprefix)
    val(label)

  output:
    path("${pprefix}_${label}_${window_size}kb.${chr}.annot.gz", emit: for_bm_part2_annot)
    path("${pprefix}_${label}_${window_size}kb.${chr}.l2.ldscore.gz", emit: for_bm_part2_ld)
    path("${pprefix}_${label}_${window_size}kb.${chr}.l2.M", emit: for_bm_part2_m)
    path("${pprefix}_${label}_${window_size}kb.${chr}.l2.M_5_50", emit: for_bm_part2_m_5_50)
    path("${pprefix}_${label}_${window_size}kb.${chr}.log", emit: for_bm_part2_log)

  script:
  """
  tail -n+2 ${gene_annot_path}|\
    awk -v OFS="\t" \
      'BEGIN{ print "ensembl_id","ensembl_strand","ensembl_bp_start","ensembl_bp_end","transcript_biotype","HGNC","ensembl_chromosome"} {if(\$6==\$5){print \$1,"-1",\$4,\$5,"protein_coding",\$2,\$3} else {print \$1,"1",\$4,\$5,"protein_coding",\$2,\$3}}' \
      > bm_gene_boundary.txt

  python ${script_dir}/bin/benchmarker/src/annotate_input_gene_list.py \
    --prio_gene_list ${prio_gene_list} \
    --curr_chromosome ${chr} \
    --gene_boundary_file bm_gene_boundary.txt \
    --output_label ${pprefix}_${label} \
    --window_size_in_kb ${window_size} \
    --baseline_file ${baseline_pathdir}/${baseline_prefix}.${chr}.annot.gz \
    --prio_gene_list_has_header ${params.benchmarker.prio_gene_list_has_header} \
    --gene_col_name ${params.benchmarker.gene_col_name}

  python ${script_dir}/bin/ldsc/ldsc.py \
    --l2 --ld-wind-cm 1 \
    --bfile ${bfile_pathdir}/${bfile_prefix}.${chr} \
    --annot ${pprefix}_${label}_${window_size}kb.${chr}.annot.gz \
    --out ${pprefix}_${label}_${window_size}kb.${chr} \
    --print-snps $print_snps
  """
}


process BENCHMARKER_PART2 {
  publishDir "${params.outdir}/5_bm_results/${pprefix}/part_2/", mode: 'copy', overwrite: true

  input:
    path(part_1_files_annot)
    path(part_1_files_ld)
    path(part_1_files_m)
    path(part_1_files_m_5_50)
    path(part_1_files_log)
    path(script_dir)
    path(munged_gwas)
    path(baseline_pathdir)
    val(baseline_prefix)
    path(gene_window_pathdir)
    val(gene_window_prefix)
    path(weights_pathdir)
    val(weights_prefix)
    path(frq_pathdir)
    val(frq_prefix)
    val(window_size)
    val(pprefix)
    val(label)

  output:
    path("${pprefix}_${label}*")

  script:
  """
  python ${script_dir}/bin/ldsc/ldsc.py \
    --overlap-annot --print-coefficients --print-delete-vals \
    --h2 ${munged_gwas} \
    --ref-ld-chr ${baseline_pathdir}/${baseline_prefix}.,${gene_window_pathdir}/${gene_window_prefix},${pprefix}_${label}_${window_size}kb. \
    --w-ld-chr ${weights_pathdir}/${weights_prefix} \
    --frqfile-chr ${frq_pathdir}/${frq_prefix} \
    --out ./${pprefix}_${label}_${window_size}kb

  python ${script_dir}/bin/benchmarker/src/add_stats_to_results_files.py \
    --trait ${pprefix} \
    --label ${label} \
    --full_file_stem ${pprefix}_${label}_${window_size}kb \
    --baseline_file_path ${baseline_pathdir}/${baseline_prefix}. \
    --input_file_directory ./ \
    --output_file_directory ./
  """
}
