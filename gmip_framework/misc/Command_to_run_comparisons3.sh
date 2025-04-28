for i in $(cat /work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/GMIP.input.csv)
do
  trait=$(echo $i|cut -d, -f1)
  magma_input=$(echo $i|cut -d, -f2)
  munged_gwas=$(echo $i|cut -d, -f3)
  gold_standard_file=$(echo $i|cut -d, -f4)
  
  for f in full_pops  expression_pops ppi_pops pathway_pops nafld
  do 
    sbatch2 ${trait}.${f} 1 1 8 MCB23017 "
    curr_path=$(pwd)
    mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_28_full/${trait}/${f}/
    cd /scratch/04179/mshabb/GMIP2/2024_09_28_full/${trait}/${f}/
    nextflow run /home1/04179/mshabb/GMIP2/workflows/1_2_3_4_5_6_magma_preprocess_pops_gather_eval.wf.nf \
      -c /home1/04179/mshabb/GMIP2/conf/featurewise/${f}.nextflow.config -profile local_conda -resume \
      --magma_input ${magma_input} \
      --munged_gwas ${munged_gwas} \
      --gold_standard_file ${gold_standard_file}
    cd ${curr_path}" "" skx
    sleep 2s
  done
  
  for f in netwas naga
  do 
    sbatch2 ${trait}.${f} 1 1 8 MCB23017 "
    curr_path=$(pwd)
    mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_28_full/${trait}/${f}/
    cd /scratch/04179/mshabb/GMIP2/2024_09_28_full/${trait}/${f}/
    nextflow run /home1/04179/mshabb/GMIP2/workflows/1_2_3_4_5_6_7_magma_preprocess_pops_naga_gather_eval.wf.nf \
      -c /home1/04179/mshabb/GMIP2/conf/featurewise/${f}.nextflow.config -profile local_conda -resume \
      --magma_input ${magma_input} \
      --munged_gwas ${munged_gwas} \
      --gold_standard_file ${gold_standard_file}
    cd ${curr_path}" "" skx
    sleep 2s
  done
done