//nextflow run ~/GMIP2/workflows/2_3_preprocess_pops_wf.nf -with-conda

params.base_refdir  = "${projectDir}/refdir"
params.script_dir   = projectDir
params.outdir       = "outdir"
params.prefix       = "out"

//preprocess related args
params.feature_mat   = null
params.magma_out     = null
params.gene_loc_file = null
params.pval_column   = "P"
params.threshold     = "5e-8"
params.chrom_column  = "CHR"
params.x_id_col      = "ENSGID"
params.y_id_col      = "GENE"
params.folds         = "3"

//pops related args
params.gene_annot_path         = "${projectDir}/bin/pops/example/data/utils/gene_annot_jun10.txt"
params.feature_mat_pathdir     = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/ppi_pops_features_munged"
params.feature_mat_prefix      = "pops_features"
params.control_features_path   = "${params.base_refdir}/2_pops_files/data_from_dropbox/ppi_pops/control.features"

//
include { PREPROCESS         }  from '../modules/preprocess.nf'

//
include { POPS as POPS_noCV  }  from '../modules/pops.nf'
include { POPS as POPS_LOCO  }  from '../modules/pops.nf'
include { POPS as POPS_kFOLD }  from '../modules/pops.nf'

//
include { GATHER_RES as GATHER_RES_POPS_noCV }  from '../modules/gather_res.nf'
include { GATHER_RES as GATHER_RES_POPS_LOCO }  from '../modules/gather_res.nf'
include { GATHER_RES as GATHER_RES_POPS_kFOLD } from '../modules/gather_res.nf'

workflow {
  PREPROCESS_otch=PREPROCESS(
    params.script_dir,
    params.pval_column,
    params.threshold,
    params.chrom_column,
    params.x_id_col,
    params.y_id_col,
    params.prefix,
    params.folds,
    params.feature_mat,
    params.magma_out,
    params.gene_loc_file
  )

  // LOCO commands
    def lineCount=0
    ml_loco_inch=PREPROCESS_otch.preprocess_chrom_csv.splitText().map { line ->
      if (lineCount++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName)
      }
    }

    POPS_LOCO_otch=POPS_LOCO(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_chrom_files.collect(), "pops_loco", ml_loco_inch
    )

    POPS_LOCO_GATHER_otch=GATHER_RES_POPS_LOCO( "pops_loco", POPS_LOCO_otch.pops_yhat.collect())

  // kFold commands
    def lineCount2=0
    ml_kfold_inch=PREPROCESS_otch.preprocess_fold_csv.splitText().map { line ->
      if (lineCount2++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName)
      }
    }

    POPS_kFOLD_otch=POPS_kFOLD(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_fold_files.collect(), "pops_kfold", ml_kfold_inch
    )

    POPS_kFOLD_GATHER_otch=GATHER_RES_POPS_kFOLD( "pops_kfold", POPS_kFOLD_otch.pops_yhat.collect())

  // noCV commands
    def lineCount3=0
    ml_nocv_inch=PREPROCESS_otch.preprocess_nocv_csv.splitText().map { line ->
      if (lineCount3++ == 0) { return null } else {
        def (fold, XtrainFileName, XtestFileName, YtrainFileName, YtestFileName, Xtest2FileName) = line.split(',')
        Xtest2FileName=Xtest2FileName.trim()
        return tuple(fold,XtrainFileName,XtestFileName,YtrainFileName,YtestFileName,Xtest2FileName)
      }
    }

    POPS_noCV_otch=POPS_noCV(
      params.script_dir, params.gene_annot_path, params.feature_mat_pathdir, params.feature_mat_prefix,
      params.control_features_path, PREPROCESS_otch.preprocess_nocv_files.collect(), "pops_nocv", ml_nocv_inch
    )

    POPS_noCV_GATHER_otch=GATHER_RES_POPS_noCV( "pops_nocv", POPS_noCV_otch.pops_yhat.collect())
}
