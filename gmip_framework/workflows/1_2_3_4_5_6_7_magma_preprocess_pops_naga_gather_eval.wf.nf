//nextflow run ~/GMIP2/workflows/1_2_3_4_5_6_magma_preprocess_pops_gather_eval.wf.nf -c ~/GMIP2/conf/NAFLD_test.nextflow.config -profile local_conda -resume

//
include { saveParams }  from '../modules/save_params.nf'
include { MAGMA }       from '../modules/magma.nf'
include { PREPROCESS }  from '../modules/preprocess.nf'

//
  include { POPS as POPS_noCV  }  from '../modules/pops.nf'
  include { POPS as POPS_LOCO  }  from '../modules/pops.nf'
  include { POPS as POPS_kFOLD }  from '../modules/pops.nf'
  
  include { NAGA as NAGA_noCV  }  from '../modules/naga.nf'
  include { NAGA as NAGA_LOCO  }  from '../modules/naga.nf'
  include { NAGA as NAGA_kFOLD }  from '../modules/naga.nf'

//
  include { GATHER_RES as GATHER_RES_POPS_noCV }  from '../modules/gather_res.nf'
  include { GATHER_RES as GATHER_RES_POPS_LOCO }  from '../modules/gather_res.nf'
  include { GATHER_RES as GATHER_RES_POPS_kFOLD } from '../modules/gather_res.nf'

  include { GATHER_RES as GATHER_RES_NAGA_noCV }  from '../modules/gather_res.nf'
  include { GATHER_RES as GATHER_RES_NAGA_LOCO }  from '../modules/gather_res.nf'
  include { GATHER_RES as GATHER_RES_NAGA_kFOLD } from '../modules/gather_res.nf'

//
  include { BENCHMARKER as BENCHMARKER_POPS_noCV_top0500  } from '../subworkflows/benchmarker.nf'
  include { BENCHMARKER as BENCHMARKER_POPS_LOCO_top0500  } from '../subworkflows/benchmarker.nf'
  include { BENCHMARKER as BENCHMARKER_POPS_kFOLD_top0500 } from '../subworkflows/benchmarker.nf'

  include { BENCHMARKER as BENCHMARKER_NAGA_noCV_top0500  } from '../subworkflows/benchmarker.nf'
  include { BENCHMARKER as BENCHMARKER_NAGA_LOCO_top0500  } from '../subworkflows/benchmarker.nf'
  include { BENCHMARKER as BENCHMARKER_NAGA_kFOLD_top0500 } from '../subworkflows/benchmarker.nf'

//
  include { EVAL as EVAL_POPS_noCV  } from '../modules/eval.nf'
  include { EVAL as EVAL_POPS_LOCO  } from '../modules/eval.nf'
  include { EVAL as EVAL_POPS_kFOLD } from '../modules/eval.nf'
  
  include { EVAL as EVAL_NAGA_noCV  } from '../modules/eval.nf'
  include { EVAL as EVAL_NAGA_LOCO  } from '../modules/eval.nf'
  include { EVAL as EVAL_NAGA_kFOLD } from '../modules/eval.nf'

workflow {
  //MAGMA
    MAGMA_otch=MAGMA(
      params.magma.executable, params.magma_input, params.magma.ref_pathdir, params.magma.bfile_prefix,
      params.gene_annot_path, params.magma.outdir, params.magma.prefix
    )

  // Preprocess
    PREPROCESS_otch=PREPROCESS(
      params.script_dir, params.preprocess.pval_column, params.preprocess.threshold, params.preprocess.chrom_column,
      params.preprocess.x_id_col, params.preprocess.y_id_col, params.prefix, params.preprocess.folds,
      params.feature_mat, MAGMA_otch.magma_out, MAGMA_otch.magma_gene_loc
    )

  // noCV commands
    def lineCount=0
    ml_nocv_inch=PREPROCESS_otch.preprocess_nocv_csv.splitText().map { line ->
      if (lineCount++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName) }
    }
    
    POPS_noCV_otch=POPS_noCV(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_nocv_files.collect(), "pops_nocv", ml_nocv_inch
    )

    NAGA_noCV_otch=NAGA_noCV(
      params.script_dir, params.feature_network, PREPROCESS_otch.preprocess_nocv_files.collect(),
      "naga_nocv", ml_nocv_inch, params.naga.ScoreColname, params.naga.IDColname
    )
    
    POPS_noCV_GATHER_otch=GATHER_RES_POPS_noCV("pops_nocv", POPS_noCV_otch.pops_yhat.collect())
    BENCHMARKER_POPS_noCV_top0500(POPS_noCV_GATHER_otch.genes_top0500, "pops_nocv", "${params.label}_top0500")
    EVAL_POPS_noCV(params.script_dir, PREPROCESS_otch.preprocess_Ytrain_full, POPS_noCV_GATHER_otch.res_all,
      "pops_nocv", params.gold_standard_file, params.eval.GENE_COL1, params.eval.SCORE_COL1, params.eval.BINARY_COL1,
      params.eval.GENE_COL2, params.eval.SCORE_COL2)

    NAGA_noCV_GATHER_otch=GATHER_RES_NAGA_noCV("naga_nocv", NAGA_noCV_otch.naga_yhat.collect())
    BENCHMARKER_NAGA_noCV_top0500(NAGA_noCV_GATHER_otch.genes_top0500, "naga_nocv", "${params.label}_top0500")
    EVAL_NAGA_noCV(params.script_dir, PREPROCESS_otch.preprocess_Ytrain_full, NAGA_noCV_GATHER_otch.res_all,
      "naga_nocv", params.gold_standard_file, params.eval.GENE_COL1, params.eval.SCORE_COL1, params.eval.BINARY_COL1,
      params.eval.GENE_COL2, params.eval.SCORE_COL2)

  // LOCO commands
    def lineCount3=0
    ml_loco_inch=PREPROCESS_otch.preprocess_chrom_csv.splitText().map { line ->
      if (lineCount3++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName) }
    }

    POPS_LOCO_otch=POPS_LOCO(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_chrom_files.collect(), "pops_loco", ml_loco_inch
    )
    
    NAGA_LOCO_otch=NAGA_LOCO(
      params.script_dir, params.feature_network, PREPROCESS_otch.preprocess_chrom_files.collect(),
      "naga_loco", ml_loco_inch, params.naga.ScoreColname, params.naga.IDColname
    )
    
    POPS_LOCO_GATHER_otch=GATHER_RES_POPS_LOCO( "pops_loco", POPS_LOCO_otch.pops_yhat.collect())
    BENCHMARKER_POPS_LOCO_top0500(POPS_LOCO_GATHER_otch.genes_top0500, "pops_loco", "${params.label}_top0500")
    EVAL_POPS_LOCO(params.script_dir, PREPROCESS_otch.preprocess_Ytrain_full, POPS_LOCO_GATHER_otch.res_all,
      "pops_loco", params.gold_standard_file, params.eval.GENE_COL1, params.eval.SCORE_COL1, params.eval.BINARY_COL1,
      params.eval.GENE_COL2, params.eval.SCORE_COL2)
  
    NAGA_LOCO_GATHER_otch=GATHER_RES_NAGA_LOCO("naga_loco", NAGA_LOCO_otch.naga_yhat.collect())
    BENCHMARKER_NAGA_LOCO_top0500(NAGA_LOCO_GATHER_otch.genes_top0500, "naga_loco", "${params.label}_top0500")
    EVAL_NAGA_LOCO(params.script_dir, PREPROCESS_otch.preprocess_Ytrain_full, NAGA_LOCO_GATHER_otch.res_all,
      "naga_loco", params.gold_standard_file, params.eval.GENE_COL1, params.eval.SCORE_COL1, params.eval.BINARY_COL1,
      params.eval.GENE_COL2, params.eval.SCORE_COL2)

  // kFold commands
    def lineCount4=0
    ml_kfold_inch=PREPROCESS_otch.preprocess_fold_csv.splitText().map { line ->
      if (lineCount4++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName) }
    }

    POPS_kFOLD_otch=POPS_kFOLD(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_fold_files.collect(), "pops_kfold", ml_kfold_inch
    )
    
    NAGA_kFOLD_otch=NAGA_kFOLD(
      params.script_dir, params.feature_network, PREPROCESS_otch.preprocess_fold_files.collect(),
      "naga_kfold", ml_kfold_inch, params.naga.ScoreColname, params.naga.IDColname
    )
    
    POPS_kFOLD_GATHER_otch=GATHER_RES_POPS_kFOLD( "pops_kfold", POPS_kFOLD_otch.pops_yhat.collect())
    BENCHMARKER_POPS_kFOLD_top0500(POPS_kFOLD_GATHER_otch.genes_top0500, "pops_kfold", "${params.label}_top0500")
    EVAL_POPS_kFOLD(params.script_dir, PREPROCESS_otch.preprocess_Ytrain_full, POPS_kFOLD_GATHER_otch.res_all,
      "pops_kfold", params.gold_standard_file, params.eval.GENE_COL1, params.eval.SCORE_COL1, params.eval.BINARY_COL1,
      params.eval.GENE_COL2, params.eval.SCORE_COL2)
    
    NAGA_kFOLD_GATHER_otch=GATHER_RES_NAGA_kFOLD("naga_kfold", NAGA_kFOLD_otch.naga_yhat.collect())
    BENCHMARKER_NAGA_kFOLD_top0500(NAGA_kFOLD_GATHER_otch.genes_top0500, "naga_kfold", "${params.label}_top0500")
    EVAL_NAGA_kFOLD(params.script_dir, PREPROCESS_otch.preprocess_Ytrain_full, NAGA_kFOLD_GATHER_otch.res_all,
      "naga_kfold", params.gold_standard_file, params.eval.GENE_COL1, params.eval.SCORE_COL1, params.eval.BINARY_COL1,
      params.eval.GENE_COL2, params.eval.SCORE_COL2)
}
