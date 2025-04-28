library(argparse)
# Parse command-line arguments
parser=ArgumentParser(description="Analyze gene lists and perform GSEA and overlap analysis")
parser$add_argument("ytestFile", help="Path to the first input file (e.g., MAGMA output)")
parser$add_argument("ypredFile", help="Path to the second input file (e.g., reprioritization scores)")
parser$add_argument("outprefix", help="Prefix for all output files")
parser$add_argument("gold_standardFile", default="none", nargs="?",
  help="Optional path to the gold_standardFile (one column gene list)")
parser$add_argument("--gene_col1", default="ENSGID", help="Column name for gene IDs in file1")
parser$add_argument("--score_col1", default="Score", help="Column name for scores in file1")
parser$add_argument("--binary_col1", default="binarized_P", help="binary column in file1")
parser$add_argument("--gene_col2", default="ENSGID", help="Column name for gene IDs in file2")
parser$add_argument("--score_col2", default="PredictedScore", help="Column name for scores in file2")
args=parser$parse_args()

library(enrichplot)
library(data.table)
library(clusterProfiler)

get_list=function(my_dt, score_col, gene_col, name_string){
  my_list=my_dt[[score_col]]
  names(my_list)=my_dt[[gene_col]]
  my_list=sort(my_list, decreasing=T)
  top_sets=c(50, 100, 250, 500, 1000)
  my_df=data.frame()
  for (n in top_sets) {
    my_df=rbind(my_df, data.frame(TERM=paste0("top_", sprintf("%04d", n), "_", name_string),
      Gene=names(my_list)[1:n]))
  }
  return(list(rankedList=my_list, term2gene=my_df))
}

INAP_eval_continuous_truth_continuous_outcome=function(y_true, y_pred){
  return(data.frame(
    MAE=MLmetrics::MAE(y_pred=y_pred, y_true=y_true),
    MAPE=MLmetrics::MAPE(y_pred=y_pred, y_true=y_true),
    MedianAE=MLmetrics::MedianAE(y_pred=y_pred, y_true=y_true),
    MedianAPE=MLmetrics::MedianAPE(y_pred=y_pred, y_true=y_true),
    MSE=MLmetrics::MSE(y_pred=y_pred, y_true=y_true),
    R2_Score=MLmetrics::R2_Score(y_pred=y_pred, y_true=y_true),
    RAE=MLmetrics::RAE(y_pred=y_pred, y_true=y_true),
    RMSE =MLmetrics::RMSE(y_pred=y_pred, y_true=y_true),
    RMSLE=MLmetrics::RMSLE(y_pred=y_pred, y_true=y_true),
    RMSPE=MLmetrics::RMSPE(y_pred=y_pred, y_true=y_true),
    PearsonCor=cor(y_pred, y_true, method = "pearson"),
    SpearmanCor=cor(y_pred, y_true, method="spearman"),
    KendallCor=cor(y_pred, y_true, method="kendall")
  ))
}

INAP_eval_binary_truth_continuous_outcome=function(y_true, y_pred){
  if(sum(y_true)!=length(y_true)){
      PRAUC=MLmetrics::PRAUC(y_pred=y_pred, y_true=y_true)
      LiftAUC=MLmetrics::LiftAUC(y_pred=y_pred, y_true=y_true)
      GainAUC=MLmetrics::GainAUC(y_pred=y_pred, y_true=y_true)
  } else { PRAUC=NaN ; LiftAUC=NaN ; GainAUC=NaN }

  return(data.frame(
    LogLoss=MLmetrics::LogLoss(y_pred=y_pred, y_true=y_true),
    Poisson_LogLoss=MLmetrics::Poisson_LogLoss(y_pred=y_pred, y_true=y_true),
    AUC=MLmetrics::AUC(y_pred=y_pred, y_true=y_true),
    Gini=MLmetrics::Gini(y_pred=y_pred, y_true=y_true),
    NormalizedGini=MLmetrics::NormalizedGini(y_pred=y_pred, y_true=y_true),
    KS_Stat=MLmetrics::KS_Stat(y_pred=y_pred, y_true=y_true),
    PRAUC=PRAUC, LiftAUC=LiftAUC, GainAUC=GainAUC
  ))
}

# Define the function to analyze the gene lists
analyze_genes=function(ytestFile, ypredFile, outprefix, gold_standardFile=NULL, gene_col1="ENSGID", score_col1="Score",
  binary_col1="binarized_P", gene_col2="ENSGID", score_col2="PredictedScore") {
  #
  ytest_dt=fread(ytestFile)
  ypred_dt=fread(ypredFile)
  merged_dt=merge(ytest_dt, ypred_dt, by.x=gene_col1, by.y=gene_col2)
  
  #
  ytest=get_list(ytest_dt, score_col1, gene_col1, "gwas_genes")
  ypred=get_list(ypred_dt, score_col2, gene_col2, "reprio_genes")

  # Calculate enrichment of predicted lists in original ranked GWAS gene list
  gsea_ytest_asRL=GSEA(ytest$rankedList, TERM2GENE=ypred$term2gene, pvalueCutoff=0.9, maxGSSize=20000, eps=1e-100)

  # Calculate enrichment of GWAS gene lists in ranked predicted gene list
  gsea_ypred_asRL=GSEA(ypred$rankedList, TERM2GENE=ytest$term2gene, pvalueCutoff=0.9, maxGSSize=20000, eps=1e-100)

  # png(paste0(outprefix, ".ytest_asRL.gsea.plot.png"), width=12, height=8, units="in", res=150)
  #   print(gseaplot2(gsea_ytest_asRL, pvalue_table=T, geneSetID=1:5,
  #     color=c("#1E90FF", "#E495A5", "#3CB371", "#DAA520", "#9932CC")))
  # dev.off()
  # 
  # png(paste0(outprefix, ".ypred_asRL.gsea.plot.png"), width=12, height=8, units="in", res=150)
  #   print(gseaplot2(gsea_ypred_asRL, pvalue_table=T, geneSetID=1:5,
  #     color=c("#1E90FF", "#E495A5", "#3CB371", "#DAA520", "#9932CC")))
  # dev.off()
  # 
  # png(paste0(outprefix, ".ytest_asRL.gsea.dotplot.png"), width=7, height=7, units="in", res=150)
  #   print(dotplot(gsea_ytest_asRL, x="NES"))
  # dev.off()
  # 
  # png(paste0(outprefix, ".ypred_asRL.gsea.dotplot.png"), width=7, height=7, units="in", res=150)
  #   print(dotplot(gsea_ypred_asRL, x="NES"))
  # dev.off()

  fwrite(gsea_ytest_asRL@result, paste0(outprefix, ".ytest_asRL.gsea.tsv"), sep="\t")
  fwrite(gsea_ypred_asRL@result, paste0(outprefix, ".ypred_asRL.gsea.tsv"), sep="\t")
  
  overlap_res_list=list()
  final_df=data.frame()
  i=1
  # Overlap analysis
  for(n in c(50, 100, 250, 500, 1000)){
    overlap_res_list[[i]]=enricher(names(ypred$rankedList)[1:n], TERM2GENE=ytest$term2gene, pvalueCutoff=1,
      universe=names(ypred$rankedList), maxGSSize=20000, qvalueCutoff=1)
    overlap_res_list[[i]]@result
    mydf=overlap_res_list[[i]]@result
    mydf$ID=paste0(mydf$ID, "__ypred_top", n, "genes")
    final_df=rbind(final_df, mydf)
    i=i+1
  }
  
  fwrite(final_df, paste0(outprefix, ".ypred_gwas.overlap.tsv"), sep="\t")

  # Optional file processing if provided
  if (gold_standardFile!="none") {
    gold_standard_genes=fread(gold_standardFile, header=FALSE)
    gs_df=data.frame(TERM="GoldStandard", Gene=gold_standard_genes[[1]])
    # Calculate enrichment of Gold Standard genes in original ranked GWAS gene list
    gsea_ytest_asRL_GS=GSEA(ytest$rankedList, TERM2GENE=gs_df, pvalueCutoff=1, maxGSSize=20000, eps=1e-100)
    fwrite(gsea_ytest_asRL_GS@result, paste0(outprefix, ".goldStandard_overlap_ytest_asRL.gsea.tsv"), sep="\t")
    # Calculate enrichment of Gold Standard genes in ranked predicted gene list
    gsea_ypred_asRL_GS=GSEA(ypred$rankedList, TERM2GENE=gs_df, pvalueCutoff=1, maxGSSize=20000, eps=1e-100)
    fwrite(gsea_ypred_asRL_GS@result, paste0(outprefix, ".goldStandard_overlap_ypred_asRL.gsea.tsv"), sep="\t")
    
    #Do overlap analysis too
    overlap_res_GS=enricher(gold_standard_genes[[1]], TERM2GENE=rbind(ytest$term2gene, ypred$term2gene),
      pvalueCutoff=1, universe=names(ytest$rankedList), maxGSSize=20000, qvalueCutoff=1)
    fwrite(overlap_res_GS@result, paste0(outprefix, ".goldStandard_ypred_gwas.overlap.tsv"), sep="\t")
  }
  
  eval_df1=INAP_eval_continuous_truth_continuous_outcome(y_true=merged_dt$Score, y_pred=merged_dt$PredictedScore)
  colnames(eval_df1)=paste0("continuous_truth.", colnames(eval_df1))
  eval_df2=INAP_eval_binary_truth_continuous_outcome(y_true=merged_dt$binarized_P, y_pred=merged_dt$PredictedScore)
  colnames(eval_df2)=paste0("binary_truth.", colnames(eval_df2))
  eval_df=rbind(t(eval_df1), t(eval_df2))
  fwrite(data.table(Measure=rownames(eval_df), Value=eval_df[,1]),
    paste0(outprefix, ".ypred_MLmeasures.tsv"), sep="\t")
  
  png(paste0(outprefix, ".goldStandard_overlap_ytest_asRL.gsea.plot.png"), width=12, height=8, units="in", res=150)
    print(gseaplot(gsea_ytest_asRL_GS,geneSetID=1))
  dev.off()
  
  png(paste0(outprefix, ".goldStandard_overlap_ypred_asRL.gsea.plot.png"), width=12, height=8, units="in", res=150)
    print(gseaplot(gsea_ypred_asRL_GS,geneSetID=1))
  dev.off()

  return("All Done")
}

# Run the function with parsed arguments
analyze_genes(
  ytestFile=args$ytestFile,
  ypredFile=args$ypredFile,
  outprefix=args$outprefix,
  gold_standardFile=args$gold_standardFile,
  gene_col1=args$gene_col1,
  score_col1=args$score_col1,
  binary_col1=args$binary_col1,
  gene_col2=args$gene_col2,
  score_col2=args$score_col2
)

ytestFile="out_Ytrain_full.tsv"
ypredFile="naga_loco.final_ytest.preds2.txt"
outprefix="naga_loco"
gold_standardFile="RA.gold_standard.genes.list.ENSGID.txt"
gene_col1="ENSGID"
score_col1="Score"
binary_col1="binarized_P"
gene_col2="ENSGID"
score_col2="PredictedScore"
# Rscript ~/GMIP2/bin/eval.R /scratch/04179/mshabb/GMIP2/new_tests/NAFLD/GMIP2/1_2_3_4_5_magma_preprocess_pops_gather_eval/outdir/2_preprocess/out_Ytrain_full.tsv.gz /scratch/04179/mshabb/GMIP2/new_tests/NAFLD/GMIP2/1_2_3_4_5_magma_preprocess_pops_gather_eval/outdir/4_gather_res/pops_nocv/pops_nocv.final_ytest.preds2.txt Test /scratch/04179/mshabb/GMIP2/data/GoldStandardGeneList/NAFLD/NAFLD.gold_standard.genes.list.ENSGID.txt