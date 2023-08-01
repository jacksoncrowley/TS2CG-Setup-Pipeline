#!/usr/bin/env nextflow

include { runTS2CG } from './modules/runts2cg.nf'
include { rotate } from './modules/rotate.nf'
include { em; em as em2; em as em3 } from './modules/em.nf'
include { eq; eq as eq2 } from './modules/eq.nf'
include { solvate } from './modules/solvate.nf'

workflow {
    // run TS2CG
    runTS2CG(params.in, params.top_header)

    // rotate TS2CG output
    rotate(runTS2CG.out.output_gro)

    // run em #1
    em(rotate.out.rotated_gro, runTS2CG.out.system_top, params.em1)

    // run em #2 
    em2(em.out.em_gro, runTS2CG.out.system_top, params.em2)

    // run eq #1
    eq(em2.out.em_gro, runTS2CG.out.system_top, params.eq1)

    // solvate
    solvate(eq.out.eq_gro, runTS2CG.out.system_top, params.water)

    // run em #3
    em3(solvate.out.solvated_gro, solvate.out.topol, params.em2)

    // run eq #2
    eq2(em3.out.em_gro, solvate.out.topol, params.eq1)
    view
}
