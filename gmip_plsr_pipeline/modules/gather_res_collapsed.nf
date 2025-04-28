process GATHER_RES_COLLAPSED {
  publishDir "${params.outdir}/gather_res/${prefix}/", mode: 'copy', overwrite: true

  input:
    val(prefix)
    val(header_string)
    path(pred2_files)
    path(loci_file)

  output:
    path("${prefix}.final_ytest.preds2.txt",               emit: res_all)
    path("${prefix}.final_ytest.preds2.all.genes.txt",     emit: genes_all)
    path("${prefix}.final_ytest.preds2.top0005.genes.txt", emit: genes_top0005)
    path("${prefix}.final_ytest.preds2.top0025.genes.txt", emit: genes_top0025)
    path("${prefix}.final_ytest.preds2.top0050.genes.txt", emit: genes_top0050)
    path("${prefix}.final_ytest.preds2.top0100.genes.txt", emit: genes_top0100)
    path("${prefix}.final_ytest.preds2.top0250.genes.txt", emit: genes_top0250)
    path("${prefix}.final_ytest.preds2.top0500.genes.txt", emit: genes_top0500)
    path("${prefix}.final_ytest.preds2.top0750.genes.txt", emit: genes_top0750)
    path("${prefix}.final_ytest.preds2.top1000.genes.txt", emit: genes_top1000)
    path("${prefix}.final_ytest.preds2.top2000.genes.txt", emit: genes_top2000)
    path("${prefix}.final_ytest.preds2.top4000.genes.txt", emit: genes_top4000)

  script:
  """
  echo ${header_string}|sed 's/,/\t/g' > ${prefix}.final_ytest.preds2.txt
  cat ${prefix}.test_chr*.preds2|sort -k2,2gr >> ${prefix}.final_ytest.preds2.txt
  tail -n+2  ${prefix}.final_ytest.preds2.txt|sort -k2,2gr|cut -f1 > ${prefix}.final_ytest.preds2.all.genes.txt
  head    -5 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0005.genes.txt
  head   -25 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0025.genes.txt
  head   -50 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0050.genes.txt
  head  -100 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0100.genes.txt
  head  -250 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0250.genes.txt
  head  -500 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0500.genes.txt
  head  -750 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top0750.genes.txt
  head -1000 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top1000.genes.txt
  head -2000 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top2000.genes.txt
  head -4000 ${prefix}.final_ytest.preds2.all.genes.txt > ${prefix}.final_ytest.preds2.top4000.genes.txt
  """
}
