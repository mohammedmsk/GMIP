trait=NAFLD
magma_input=/work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/2_magma_input/NAFLD/miao_etal_hgg/NAFLD_imp_bgen.magma_input.tsv
munged_gwas=/work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/3_munged_gwas/NAFLD/miao_etal_hgg/NAFLD_imp_bgen.forLDSC.sumstats.gz
gold_standard_file=/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/data/GoldStandardGeneList/NAFLD/NAFLD.gold_standard.genes.list.ENSGID.txt

for f in full_pops  expression_pops ppi_pops pathway_pops nafld
do 
  sbatch2 ${trait}.${f} 1 1 8 MCB23017 "
  curr_path=$(pwd)
  mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  cd /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
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
  mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  cd /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  nextflow run /home1/04179/mshabb/GMIP2/workflows/1_2_3_4_5_6_7_magma_preprocess_pops_naga_gather_eval.wf.nf \
    -c /home1/04179/mshabb/GMIP2/conf/featurewise/${f}.nextflow.config -profile local_conda -resume \
    --magma_input ${magma_input} \
    --munged_gwas ${munged_gwas} \
    --gold_standard_file ${gold_standard_file}
  cd ${curr_path}" "" skx
  sleep 2s
done

trait=RA
magma_input=/work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/2_magma_input/pops_paper/RA_GWASmeta_European_v2.magma_input.txt
munged_gwas=/work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/3_munged_gwas/pops_paper/RA_GWASmeta_European_v2.forLDSC.sumstats.gz
gold_standard_file=/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/data/GoldStandardGeneList/RA/RA.gold_standard.genes.list.ENSGID.txt

for f in full_pops  expression_pops ppi_pops pathway_pops nafld
do 
  sbatch2 ${trait}.${f} 1 1 8 MCB23017 "
  curr_path=$(pwd)
  mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  cd /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
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
  mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  cd /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  nextflow run /home1/04179/mshabb/GMIP2/workflows/1_2_3_4_5_6_7_magma_preprocess_pops_naga_gather_eval.wf.nf \
    -c /home1/04179/mshabb/GMIP2/conf/featurewise/${f}.nextflow.config -profile local_conda -resume \
    --magma_input ${magma_input} \
    --munged_gwas ${munged_gwas} \
    --gold_standard_file ${gold_standard_file}
  cd ${curr_path}" "" skx
  sleep 2s
done

trait=SCZ
magma_input=/work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/2_magma_input/benchmarker_paper/daner_PGC_SCZ52_0513a.hq2.magma_input.tsv
munged_gwas=/work2/04179/mshabb/common/GNAP_v4.0_2024_03_19/data/3_munged_gwas/benchmarker_paper/daner_PGC_SCZ52_0513a.hq2.forLDSC.sumstats.gz
gold_standard_file=/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/data/GoldStandardGeneList/SCZ/SCZ.gold_standard.genes.list.ENSGID.txt

for f in full_pops  expression_pops ppi_pops pathway_pops nafld
do 
  sbatch2 ${trait}.${f} 1 1 8 MCB23017 "
  curr_path=$(pwd)
  mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  cd /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
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
  mkdir -p /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  cd /scratch/04179/mshabb/GMIP2/2024_09_27/${trait}/${f}/
  nextflow run /home1/04179/mshabb/GMIP2/workflows/1_2_3_4_5_6_7_magma_preprocess_pops_naga_gather_eval.wf.nf \
    -c /home1/04179/mshabb/GMIP2/conf/featurewise/${f}.nextflow.config -profile local_conda -resume \
    --magma_input ${magma_input} \
    --munged_gwas ${munged_gwas} \
    --gold_standard_file ${gold_standard_file}
  cd ${curr_path}" "" skx
  sleep 2s
done

