// Include the benchmarker subworkflows
include { BENCHMARKER as BENCHMARKER } from '../subworkflows/benchmarker.2.nf'
//
workflow {
  BENCHMARKER(params.gene_list, params.prefix, params.label)
}


