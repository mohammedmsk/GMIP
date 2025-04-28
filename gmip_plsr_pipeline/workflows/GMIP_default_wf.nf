//
include { saveParams }  from '../modules/save_params.nf'

//
include { MAGMA } from '../modules/magma.nf'
include { POPS }  from '../modules/pops.nf'
include { GMIP }  from '../modules/gmip.nf'

//
include { COLLAPSE_GENES_TO_LOCI } from "../modules/collapse_genes_to_loci.nf"

//
include { GATHER_RES as GATHER_RES_POPS } from '../modules/gather_res.nf'
include { GATHER_RES as GATHER_RES_GMIP } from '../modules/gather_res.nf'

//
include { BENCHMARKER as BENCHMARKER_POPS_top0005 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0005 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top0025 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0025 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top0050 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0050 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top0100 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0100 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top0250 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0250 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top0500 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0500 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top0750 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top0750 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top1000 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top1000 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top2000 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top2000 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_POPS_top4000 } from '../subworkflows/benchmarker.nf'
include { BENCHMARKER as BENCHMARKER_GMIP_top4000 } from '../subworkflows/benchmarker.nf'

//
workflow {

  MAGMA_otch=MAGMA(
    params.magma.executable,
    params.magma_input,
    params.magma.ref_pathdir,
    params.magma.bfile_prefix,
    params.gene_annot_path,
    params.magma.outdir,
    params.magma.prefix
  )

  // Define a channel with all autosomal chromosome numbers
  chromosome_ch=Channel.from(1..22)

  POPS_inch=chromosome_ch.map { chr ->
    tuple(
      params.script_dir,
      params.gene_annot_path,
      params.feature_mat_pathdir,
      params.feature_mat_prefix,
      params.control_features_path,
      chr)
  }

  if (params.collapse_loci) {
    // First collapsing neighbouring genes to loci
    COLLAPSE_GENES_TO_LOCI_otch=COLLAPSE_GENES_TO_LOCI(
      MAGMA_otch.magma_out,
      params.collapse_genes_to_loci.id,
      params.collapse_genes_to_loci.chr,
      params.collapse_genes_to_loci.start,
      params.collapse_genes_to_loci.end,
      params.collapse_genes_to_loci.score,
      params.collapse_genes_to_loci.wd_sz
    )
    // Then running pops on new collapsed loci representative gene
    POPS_otch=POPS(
      POPS_inch,
      COLLAPSE_GENES_TO_LOCI_otch.magma_out,
      COLLAPSE_GENES_TO_LOCI_otch.magma_raw)
  } else {
    // No collapsing, simply running pops on magma out.
    POPS_otch=POPS(
      POPS_inch,
      MAGMA_otch.magma_out,
      MAGMA_otch.magma_raw)
  }

  GMIP_otch=GMIP(
    params.script_dir,
    params.gene_annot_path,
    POPS_otch.pops_traindata,
    POPS_otch.pops_matdata,
    POPS_otch.pops_marginals,
    params.gmip.prefix,
    params.gmip_method,
    params.gmip.noFeatureClustering
  )

  if (params.collapse_loci) {
    POPS_GATHER_otch=GATHER_RES_POPS(
      params.pops.prefix,
      "ENSGID,PoPS_Score,Y,feature_selection_gene,training_gene",
      POPS_otch.pops_preds2.collect())

    GMIP_GATHER_otch=GATHER_RES_GMIP(
      params.gmip.prefix,
      params.gmip.pred_header_string,
      GMIP_otch.gmip_preds2.collect()
    )
  } else {
    POPS_GATHER_otch=GATHER_RES_POPS(
      params.pops.prefix,
      params.pops.pred_header_string,
      POPS_otch.pops_preds2.collect())
    GMIP_GATHER_otch=GATHER_RES_GMIP(
      params.gmip.prefix,
      params.gmip.pred_header_string,
      GMIP_otch.gmip_preds2.collect()
    )
  }

  //
  BENCHMARKER_POPS_top0005(POPS_GATHER_otch.genes_top0005, params.pops.prefix, "${params.label}_top0005")
  BENCHMARKER_POPS_top0025(POPS_GATHER_otch.genes_top0025, params.pops.prefix, "${params.label}_top0025")
  BENCHMARKER_POPS_top0050(POPS_GATHER_otch.genes_top0050, params.pops.prefix, "${params.label}_top0050")
  BENCHMARKER_POPS_top0100(POPS_GATHER_otch.genes_top0100, params.pops.prefix, "${params.label}_top0100")
  BENCHMARKER_POPS_top0250(POPS_GATHER_otch.genes_top0250, params.pops.prefix, "${params.label}_top0250")
  BENCHMARKER_POPS_top0500(POPS_GATHER_otch.genes_top0500, params.pops.prefix, "${params.label}_top0500")
  BENCHMARKER_POPS_top0750(POPS_GATHER_otch.genes_top0750, params.pops.prefix, "${params.label}_top0750")
  BENCHMARKER_POPS_top1000(POPS_GATHER_otch.genes_top1000, params.pops.prefix, "${params.label}_top1000")
  BENCHMARKER_POPS_top2000(POPS_GATHER_otch.genes_top2000, params.pops.prefix, "${params.label}_top2000")
  BENCHMARKER_POPS_top4000(POPS_GATHER_otch.genes_top4000, params.pops.prefix, "${params.label}_top4000")

  BENCHMARKER_GMIP_top0005(GMIP_GATHER_otch.genes_top0005, params.gmip.prefix, "${params.label}_top0005")
  BENCHMARKER_GMIP_top0025(GMIP_GATHER_otch.genes_top0025, params.gmip.prefix, "${params.label}_top0025")
  BENCHMARKER_GMIP_top0050(GMIP_GATHER_otch.genes_top0050, params.gmip.prefix, "${params.label}_top0050")
  BENCHMARKER_GMIP_top0100(GMIP_GATHER_otch.genes_top0100, params.gmip.prefix, "${params.label}_top0100")
  BENCHMARKER_GMIP_top0250(GMIP_GATHER_otch.genes_top0250, params.gmip.prefix, "${params.label}_top0250")
  BENCHMARKER_GMIP_top0500(GMIP_GATHER_otch.genes_top0500, params.gmip.prefix, "${params.label}_top0500")
  BENCHMARKER_GMIP_top0750(GMIP_GATHER_otch.genes_top0750, params.gmip.prefix, "${params.label}_top0750")
  BENCHMARKER_GMIP_top1000(GMIP_GATHER_otch.genes_top1000, params.gmip.prefix, "${params.label}_top1000")
  BENCHMARKER_GMIP_top2000(GMIP_GATHER_otch.genes_top2000, params.gmip.prefix, "${params.label}_top2000")
  BENCHMARKER_GMIP_top4000(GMIP_GATHER_otch.genes_top4000, params.gmip.prefix, "${params.label}_top4000")
  //
  saveParams()
}
