#!/usr/bin/env nextflow

include { runTS2CG } from './modules/runts2cg.nf'
include { rotate } from './modules/rotate.nf'
include { em; em as em2 } from './modules/em.nf'

workflow {
    // run TS2CG
    runTS2CG(params.in, params.top_header)

    // rotate TS2CG output
    rotate(runTS2CG.out.output_gro)

    // run em #1
    em(rotate.out.rotated_gro, runTS2CG.out.system_top, params.em1)

    // run em.2 
    em2(em.out.em_gro, runTS2CG.out.system_top, params.em2)
    view
}
