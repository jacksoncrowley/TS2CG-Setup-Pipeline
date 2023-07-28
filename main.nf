#!/usr/bin/env nextflow

include { runTS2CG } from './modules/runts2cg.nf'
include { rotate } from './modules/rotate.nf'

workflow {

    // run TS2CG
    Channel.fromPath(params.in) | runTS2CG

    // rotate TS2CG output
    rotate(runTS2CG.out.output_gro)

    // run em
    // em(rotate.out.rotated_gro, runTS2CG.out.output_top, params.em1)
    view
}
