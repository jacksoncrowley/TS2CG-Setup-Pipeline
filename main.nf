#!/usr/bin/env nextflow

// paths to required input files
params.in = "$baseDir/generate.str"
params.em1 = "$baseDir/mdp/em1.mdp"
params.em2 = "$baseDir/mdp/em2.mdp"
params.eq1 = "$baseDir/mdp/eq1.mdp"
params.top_header = "$baseDir/top/header.txt"
params.rotate = [0, 0, 0]

// number of processes to be used by energy minimization and equilibration
params.cores = 16

// import modules
include { runTS2CG } from './modules/runts2cg.nf'
include { rotate } from './modules/rotate.nf'
include { em; em as em2; em as em3 } from './modules/em.nf'
include { eq; eq as eq2 } from './modules/eq.nf'
include { solvate } from './modules/solvate.nf'


workflow {
    // run TS2CG
    runTS2CG(params.in, params.top_header)

    // Rotate TS2CG output, if necessary
    if ( params.rotate != [0, 0, 0] ) {
    // ensure the provided list is 3 elements of integers between 0 and 360
    assert params.rotate.size() == 3 && params.rotate.every { it in 0..360 }
    // convert the list to a string, then pass the string to the rotate process
    rotate_str = params.rotate.join(" ")
    rotate(runTS2CG.out.output_gro, rotate_str)
    // run em #1 on the rotated file
    em(rotate.out.rotated_gro, runTS2CG.out.system_top, params.em1)
    } 
    else {
    // run em #1 on the TS2CG output
    em(runTS2CG.out.output_gro, runTS2CG.out.system_top, params.em1)  
    }
 
    // run em #2 
    em2(em.out.em_gro, runTS2CG.out.system_top, params.em2)

    // run eq #1
    eq(em2.out.em_gro, runTS2CG.out.system_top, params.eq1)

    // solvate
    solvate(eq.out.eq_gro, runTS2CG.out.system_top)

    // run em #3
    em3(solvate.out.solvated_gro, solvate.out.topol, params.em2)

    // run eq #2
    eq2(em3.out.em_gro, solvate.out.topol, params.eq1)
    view
}
