process MAGMA {
  publishDir "${outdir}", mode: 'copy', overwrite: true

  input:
    path(magma_executable)
    path(input_file)
    path(ref_pathdir)
    val(bfile_prefix)
    path(gene_annot_path)
    val(outdir)
    val(prefix)


  output:
    path("${prefix}.genes.out",  emit: magma_out)
    path("${prefix}.genes.raw",  emit: magma_raw)
    path("gene_loc.txt")
    path("magma_0kb.genes.annot")

  script:
  """
  cut -f1,3,4,5 ${gene_annot_path} > gene_loc.txt
  
  ./${magma_executable} \
    --annotate window=0 \
    --snp-loc ${ref_pathdir}/${bfile_prefix}.bim \
    --gene-loc gene_loc.txt \
    --out magma_0kb
    
  ./${magma_executable} \
    --bfile ${ref_pathdir}/${bfile_prefix} \
    --pval ${input_file} ncol=N \
    --gene-annot magma_0kb.genes.annot \
    --out ${prefix} \
    --gene-model snp-wise=mean
  """
}
