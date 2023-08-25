#!/usr/bin/env nextflow

// paths to required input files
params.pcg          = "$baseDir/PCG"            // PCG from TS2CG
params.in           = "$baseDir/generate.str"   // input .str file
params.outDir       = "$baseDir/results"
params.em1          = "$baseDir/mdp/em1.mdp"
params.em2          = "$baseDir/mdp/em2.mdp"
params.eq1          = "$baseDir/mdp/eq1.mdp"
params.top_header   = "$baseDir/top/header.txt" // file listing .itps in order


// take rotate parameters as string, convert to list
params.rotate = "0,0,0"
params.rotatelist = params.rotate.split(',').collect { it.toInteger() }

// number of processes to be used by energy minimization and equilibration
params.cores = 16

// import modules
include { RUNTS2CG } from './modules/runts2cg.nf'
include { ROTATE } from './modules/rotate.nf'
include { EM; EM as EM2; EM as EM3 } from './modules/em.nf'
include { EQ; EQ as EQ2 } from './modules/eq.nf'
include { SOLVATE } from './modules/solvate.nf'
include { MAKEINDEX } from './modules/make_index.nf'


workflow {
    // run TS2CG
    RUNTS2CG(params.pcg, params.in, params.top_header)

    // Rotate TS2CG output, if necessary
    if ( params.rotatelist != [0, 0, 0] ) {
        // ensure the provided list is 3 elements of integers between 0 and 360
        assert params.rotatelist.size() == 3 && params.rotatelist.every { it in 0..360 }
        // convert the list to a string, then pass the string to the rotate process
        rotate_str = params.rotatelist.join(" ")
        ROTATE(RUNTS2CG.out.output_gro, rotate_str)
        // run EM #1 on the rotated file
        EM(ROTATE.out.rotated_gro, RUNTS2CG.out.system_top, params.em1)
    } 
    else {
        // run EM #1 on the TS2CG output
        EM(RUNTS2CG.out.output_gro, RUNTS2CG.out.system_top, params.em1)  
    }
 
    // run EM #2 
    EM2(EM.out.em_gro, RUNTS2CG.out.system_top, params.em2)

    // run EQ #1
    EQ(EM2.out.em_gro, RUNTS2CG.out.system_top, params.eq1)

    // solvate
    SOLVATE(EQ.out.eq_gro, RUNTS2CG.out.system_top)

    // create index of System, Solvent, Solute
    MAKEINDEX(SOLVATE.out.solvated_gro)

    // run EM #3
    EM3(SOLVATE.out.solvated_gro, SOLVATE.out.topol, params.em2)

    // run eq #2
    EQ2(EM3.out.em_gro, SOLVATE.out.topol, params.eq1)
    
 
}
