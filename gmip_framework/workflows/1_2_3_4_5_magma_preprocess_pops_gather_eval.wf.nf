//nextflow run ~/GMIP2/workflows/1_2_3_4_5_magma_preprocess_pops_gather_eval.wf.nf -c ~/GMIP2/conf/NAFLD_test.nextflow.config -profile local_conda -resume

//
include { saveParams }  from '../modules/save_params.nf'
include { MAGMA }       from '../modules/magma.nf'
include { PREPROCESS }  from '../modules/preprocess.nf'

//
include { POPS as POPS_noCV  }  from '../modules/pops.nf'
include { POPS as POPS_LOCO  }  from '../modules/pops.nf'
include { POPS as POPS_kFOLD }  from '../modules/pops.nf'

//
include { GATHER_RES as GATHER_RES_POPS_noCV }  from '../modules/gather_res.nf'
include { GATHER_RES as GATHER_RES_POPS_LOCO }  from '../modules/gather_res.nf'
include { GATHER_RES as GATHER_RES_POPS_kFOLD } from '../modules/gather_res.nf'

//
include { BENCHMARKER as BENCHMARKER_POPS_noCV_top0500  } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_LOCO_top0500  } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_kFOLD_top0500 } from '../subworkflows/benchmarker.nf'

workflow {
  //MAGMA
    MAGMA_otch=MAGMA(
      params.magma.executable, params.magma_input, params.magma.ref_pathdir, params.magma.bfile_prefix,
      params.gene_annot_path, params.magma.outdir, params.magma.prefix
    )

  // Preprocess
    PREPROCESS_otch=PREPROCESS(
      params.script_dir, params.pval_column, params.threshold, params.chrom_column, params.x_id_col, params.y_id_col,
      params.prefix, params.folds, params.feature_mat, MAGMA_otch.magma_out, MAGMA_otch.magma_gene_loc
    )

  // LOCO commands
    def lineCount=0
    ml_loco_inch=PREPROCESS_otch.preprocess_chrom_csv.splitText().map { line ->
      if (lineCount++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName) }
    }

    POPS_LOCO_otch=POPS_LOCO(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_chrom_files.collect(), "pops_loco", ml_loco_inch
    )
    POPS_LOCO_GATHER_otch=GATHER_RES_POPS_LOCO( "pops_loco", POPS_LOCO_otch.pops_yhat.collect())
    BENCHMARKER_POPS_LOCO_top0500(POPS_LOCO_GATHER_otch.genes_top0500, "pops_loco", "${params.label}_top0500")

  // kFold commands
    def lineCount2=0
    ml_kfold_inch=PREPROCESS_otch.preprocess_fold_csv.splitText().map { line ->
      if (lineCount2++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName) }
    }

    POPS_kFOLD_otch=POPS_kFOLD(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_fold_files.collect(), "pops_kfold", ml_kfold_inch
    )

    POPS_kFOLD_GATHER_otch=GATHER_RES_POPS_kFOLD( "pops_kfold", POPS_kFOLD_otch.pops_yhat.collect())
    BENCHMARKER_POPS_kFOLD_top0500(POPS_kFOLD_GATHER_otch.genes_top0500, "pops_kfold", "${params.label}_top0500")

  // noCV commands
    def lineCount3=0
    ml_nocv_inch=PREPROCESS_otch.preprocess_nocv_csv.splitText().map { line ->
      if (lineCount3++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName) }
    }

    POPS_noCV_otch=POPS_noCV(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_nocv_files.collect(), "pops_nocv", ml_nocv_inch
    )

    POPS_noCV_GATHER_otch=GATHER_RES_POPS_noCV( "pops_nocv", POPS_noCV_otch.pops_yhat.collect())
    BENCHMARKER_POPS_noCV_top0500(POPS_noCV_GATHER_otch.genes_top0500, "pops_nocv", "${params.label}_top0500")

}
