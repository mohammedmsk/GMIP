## GWAS

We collected GWAS summary statistics for 4 traits trying to choose one from each method which will be compared.

1.  AMD GWAS from NetWAS: (<https://csg.sph.umich.edu/abecasis/public/amdgene2012/>)

```
export GMIP2_base=/scratch/04179/mshabb/GMIP2/
mkdir -p ${GMIP2_base}/data/GWAS/AMD/
cd ${GMIP2_base}/data/GWAS/AMD/

wget -c https://csg.sph.umich.edu/abecasis/public/amdgene2012/AMDGene2013_Advanced_v_Controls.txt.gz

zcat AMDGene2013_Advanced_v_Controls.txt.gz|\
  awk 'NR>1{print $1"\t"$6"\t"$4+$5}'|\
  sed '1s/^/SNP\tP\tN\n/' \
  > 2_magma_input.AMD.tsv

mkdir -p ${GMIP2_base}/data/MAGMA/AMD/
cd ${GMIP2_base}/data/MAGMA/AMD/

nextflow run ~/GMIP/workflows/GMIP_magma_only_wf.nf \
  -c ~/GMIP/conf/GMIP_magma_only.nextflow.config \
  --magma_input ${GMIP2_base}/data/GWAS/AMD/2_magma_input.AMD.tsv \
  --prefix AMD --outdir outdir -profile local_conda -resume

mkdir -p ${GMIP2_base}/data/GoldStandardGeneList/AMD
cd ${GMIP2_base}/data/GoldStandardGeneList/AMD
ls OMIM-Gene-Map-Retrieval_advanced_age-related_macular_degeneration.tsv
tail -n+6 OMIM-Gene-Map-Retrieval_advanced_age-related_macular_degeneration.tsv |cut -f6|sed 's/ //g'|tr ',' '\n'|sort|uniq > AMD.gold_standard.genes.list.txt

#Above list comes from OMIM database, gene map table for keyword search advanced age-related macular degeneration
#Above file was also manually cleaned to remove some unwanted enteries which were not gene names.

```

2.  Schizophrenia from NAGA

```
mkdir -p ${GMIP2_base}/data/GWAS/SCZ/
cd ${GMIP2_base}/data/GWAS/SCZ/

#This is the link to schizophrenia paper https://www.nature.com/articles/ng.940 Schizoprhenia paper
#This from NAGA paper 9,394 cases and 12,462 controls, so N = 21856
wget -c http://nbgwas.ucsd.edu/nagadata/schizophrenia.txt
awk 'NR>1{print $1"\t"$8"\t21856"}' schizophrenia.txt|sed '1s/^/SNP\tP\tN\n/' > 2_magma_input.SCZ.tsv

mkdir -p ${GMIP2_base}/data/MAGMA/SCZ/
cd ${GMIP2_base}/data/MAGMA/SCZ/

nextflow run ~/GMIP/workflows/GMIP_magma_only_wf.nf \
  -c ~/GMIP/conf/GMIP_magma_only.nextflow.config \
  --magma_input ${GMIP2_base}/data/GWAS/SCZ/2_magma_input.SCZ.tsv \
  --prefix SCZ --outdir outdir -profile local_conda -resume

mkdir -p ${GMIP2_base}/data/GoldStandardGeneList/SCZ
cd ${GMIP2_base}/data/GoldStandardGeneList/SCZ
ls SCZ.gold_standard.genes.list.txt
#Above gene list has been extracted from http://www.szgene.org/ as indicated in Naga paper.
#Above Only has genes corresponding to chr 1-22.
#Original paper describes to have 1147 genes but these are 816.

```

3.  NAFDL from GMIP

Website to go to get the files: https://doi.org/10.25346/S6/ENLMWL
Original GWAS paper: https://pubmed.ncbi.nlm.nih.gov/35047847/


```
mkdir -p ${GMIP2_base}/data/GWAS/NAFLD/
cd ${GMIP2_base}/data/GWAS/NAFLD

#We get the below file which has GWAS statistics. And now we convert it into format that MAGMA needs.
#We choose P_BOLT_LMM_INF as p-value as the manuscript also uses it
#N=137048 comes from the manuscript, where 28,396 NAFLD cases and 108,652 healthy individuals
zcat NAFLD_imp_bgen.stats.gz|awk 'NR>1{print $1"\t"$14"\t137048"}'|sed '1s/^/SNP\tP\tN\n/' > 2_magma_input.NAFLD.tsv

mkdir -p ${GMIP2_base}/data/MAGMA/NAFLD/
cd ${GMIP2_base}/data/MAGMA/NAFLD

nextflow run ~/GMIP/workflows/GMIP_magma_only_wf.nf \
  -c ~/GMIP/conf/GMIP_magma_only.nextflow.config \
  --magma_input ${GMIP2_base}/data/GWAS/NAFLD/2_magma_input.NAFLD.tsv \
  --prefix NAFLD --outdir outdir -profile local_conda -resume

mkdir -p ${GMIP2_base}/data/GoldStandardGeneList/NAFLD
cd ${GMIP2_base}/data/GoldStandardGeneList/NAFLD/
wc -l WP4396-datanodes.tsv
cut -f7 WP4396-datanodes.tsv|grep -v '^$'|sed 's/hgnc.symbol://'|tail -n+2|sort|uniq > NAFLD.gold_standard.genes.list.txt

#Above list comes from WikiPathways database, keyword searched:Nonalcoholic fatty liver disease
#Then selecting https://www.wikipathways.org/pathways/WP4396.html

```

4.  RA from PoPS

```
mkdir -p ${GMIP2_base}/data/GWAS/RA/
cd ${GMIP2_base}/data/GWAS/RA

wget -c https://grasp.nhlbi.nih.gov/downloads/ResultsOctober2016/Okada/RA_GWASmeta_European_v2.txt.gz
zcat RA_GWASmeta_European_v2.txt.gz|grep '^rs'|awk '{print $1"\t"$9"\t57284"}'|sed '1s/^/SNP\tP\tN\n/' > 2_magma_input.RA.tsv

mkdir -p ${GMIP2_base}/data/MAGMA/RA/
cd ${GMIP2_base}/data/MAGMA/RA

nextflow run ~/GMIP/workflows/GMIP_magma_only_wf.nf \
  -c ~/GMIP/conf/GMIP_magma_only.nextflow.config \
  --magma_input ${GMIP2_base}/data/GWAS/RA/2_magma_input.RA.tsv \
  --prefix RA --outdir outdir -profile local_conda -resume

mkdir -p ${GMIP2_base}/data/GoldStandardGeneList/RA
cd ${GMIP2_base}/data/GoldStandardGeneList/RA/
ls extracted_gene_symbols.txt
#Above list comes from https://maayanlab.cloud/Harmonizome/gene_set/Rheumatoid+Arthritis/DisGeNET+Gene-Disease+Associations
#This is disgenet from harmonizome database
#Direct download from latest disgenet website requires subscription.

```
