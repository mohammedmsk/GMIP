#cd /scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files
#mkdir -p features_munged/
#/work2/04179/mshabb/ls6/miniforge3/envs/GNAP/bin/python ~/GMIP/bin/pops/munge_feature_directory.py --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt --feature_dir data_from_dropbox/feature_dir/ --save_prefix features_munged_expression/pops_features --max_cols 100

#cd /scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox
#mkdir feature_dir_expression
#R

pops_features=data.table::fread("/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir/PoPS.features.txt.gz")
dim(pops_features)
pops_features_meta=data.table::fread("/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/gene_features_metadata.txt")

#Get expression only features
select_cols=colnames(pops_features)[colnames(pops_features) %in% c('ENSGID', pops_features_meta[pops_features_meta$V4=="Expression",]$V1)]
pops_features_expression=pops_features[,..select_cols]
dim(pops_features_expression)
pops_features_expression[1:3, 1:3]
dir.create("/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir_expression/", recursive=T)
data.table::fwrite(pops_features_expression, "/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir_expression/Expression_only.PoPS.features.txt.gz",
 sep="\t")

#Get Pathway only features
select_cols=colnames(pops_features)[colnames(pops_features) %in% c('ENSGID', pops_features_meta[pops_features_meta$V4=="Pathway",]$V1)]
pops_features_pathway=pops_features[,..select_cols]
dim(pops_features_pathway)
pops_features_pathway[1:3, 1:3]
dir.create("/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir_pathway/", recursive=T)
data.table::fwrite(pops_features_pathway, "/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir_pathway/pathway_only.PoPS.features.txt.gz", sep="\t")

#Get Pathway only features
select_cols=colnames(pops_features)[colnames(pops_features) %in% c('ENSGID', pops_features_meta[pops_features_meta$V4=="PPI",]$V1)]
pops_features_ppi=pops_features[,..select_cols]
dim(pops_features_ppi)
pops_features_ppi[1:3, 1:3]
dir.create("/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir_ppi/", recursive=T)
data.table::fwrite(pops_features_ppi, "/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files/data_from_dropbox/feature_dir_ppi/ppi_only.PoPS.features.txt.gz", sep="\t")

#cd /scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_03_19/refdir/2_pops_files
#mkdir -p features_munged_expression/
#/work2/04179/mshabb/ls6/miniforge3/envs/GNAP/bin/python ~/GMIP/bin/pops/munge_feature_directory.py --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt --feature_dir data_from_dropbox/feature_dir_expression/ --save_prefix features_munged_expression/pops_features --max_cols 100

#mkdir -p features_munged_pathway/
#/work2/04179/mshabb/ls6/miniforge3/envs/GNAP/bin/python ~/GMIP/bin/pops/munge_feature_directory.py --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt --feature_dir data_from_dropbox/feature_dir_pathway/ --save_prefix features_munged_pathway/pops_features --max_cols 100

#mkdir -p features_munged_ppi/
#/work2/04179/mshabb/ls6/miniforge3/envs/GNAP/bin/python ~/GMIP/bin/pops/munge_feature_directory.py --gene_annot_path ~/GMIP/bin/pops/example/data/utils/gene_annot_jun10.txt --feature_dir data_from_dropbox/feature_dir_ppi/ --save_prefix features_munged_ppi/pops_features --max_cols 100

#Commands for munging stats
for i in blood_WHITE_COUNT.sumstats.gz;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/') ;done
for i in LDL_with_Effect.tbl;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/') --ignore GC.Zscore;done
for i in HDL_with_Effect.tbl;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/') --ignore GC.Zscore;done
for i in TC_with_Effect.tbl;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/') --ignore GC.Zscore;done
for i in TG_with_Effect.tbl;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/') --ignore GC.Zscore;done
for i in body_BMIz.sumstats.gz body_WHRadjBMIz.sumstats.gz bp_DIASTOLICadjMEDz.sumstats.gz bp_SYSTOLICadjMEDz.sumstats.gz cov_EDU_YEARS.sumstats.gz cov_SMOKING_STATUS.sumstats.gz disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz disease_T2D.sumstats.gz pigment_SKIN.sumstats.gz repro_MENARCHE_AGE.sumstats.gz repro_MENOPAUSE_AGE.sumstats.gz;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/');done
for i in height____GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/');done
for i in daner_PGC_SCZ52_0513a.hq2.gz;do ${ldsc_python} ~/GNAP/ldsc/munge_sumstats.py --sumstats $i --out ../2_munge_sumstats/$(echo $i|sed 's/.sumstats.gz/.forLDSC/') --daner;done

#For NAFLD
/work2/04179/mshabb/ls6/miniforge3/envs/GNAP_ldsc/bin/python ~/GMIP/bin/ldsc/munge_sumstats.py --sumstats NAFLD_imp_bgen.stats.gz --out ../../../3_munged_gwas/NAFLD/miao_etal_hgg/NAFLD_imp_bgen.forLDSC --N-cas 28396 --N-con 108652 --p P_BOLT_LMM_INF --signed-sumstats BETA,0 --a1 ALLELE0 --a2 ALLELE1 --frq A1FREQ
zcat NAFLD_imp_bgen.stats.gz|cut -f1,14|awk 'NR>1{print $0"\t137048"}'|sed '1s/^/SNP\tP\tN\n/' > ../../../2_magma_input/NAFLD/miao_etal_hgg/NAFLD_imp_bgen.magma_input.tsv





























