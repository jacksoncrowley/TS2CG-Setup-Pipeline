#!/usr/bin/env nextflow

// paths to required input files
params.pcg          = "$baseDir/PCG"            // PCG from TS2CG
params.input        = "$baseDir/generate.str"   // input .str file
params.outDir       = "$baseDir/results"
params.em1          = "$baseDir/mdp/em1.mdp"
params.em2          = "$baseDir/mdp/em2.mdp"
params.eq1          = "$baseDir/mdp/eq1.mdp"
params.top_header   = "$baseDir/top/header.txt" // file listing .itps in order

// pore formation parameters
params.createporego = "$baseDir/scripts/create_pore.go"
params.poreradius   = 0
params.poreaxis     = "z"
//params.porecenter = "0,0"
//params.poreopen   = false

params.buildtop_python = "$baseDir/scripts/build_top.py"


// take rotate parameters as string, convert to list
params.rotate = "0,0,0"
params.rotatelist = params.rotate.split(',').collect { it.toInteger() }


// number of processes to be used by energy minimization and equilibration
params.cores = 16

// import modules
include { RUNTS2CG } from './modules/runts2cg.nf'
// include { ROTATE } from './modules/rotate.nf'
// include { CREATEPORE} from './modules/create_pore.nf'
// include { BUILDTOP } from './modules/build_top.nf'
include { EM; EM as EM2; EM as EM3 } from './modules/em.nf'
include { EQ; EQ as EQ2 } from './modules/eq.nf'
include { SOLVATE } from './modules/solvate.nf'
include { MAKEINDEX } from './modules/make_index.nf'


workflow {
    // run TS2CG
    RUNTS2CG(params.pcg, params.input, params.top_header)
    working_gro = RUNTS2CG.out.output_gro
    working_top = RUNTS2CG.out.system_top

    // // Steps are not yet fully developed
    // // Rotate TS2CG output, if necessary
    // if ( params.rotatelist != [0, 0, 0] ) {
    //     // ensure the provided list is 3 elements of integers between 0 and 360
    //     assert params.rotatelist.size() == 3 && params.rotatelist.every { it in 0..360 }
    //     // convert the list to a string, then pass the string to the rotate process
    //     rotate_str = params.rotatelist.join(" ")
    //     ROTATE(working_gro, rotate_str)
    //     working_gro = ROTATE.out.output_gro   
    // } 
    

    // // add pore for water EQ if needed (Vesicles)
    // if ( params.poreradius != 0 ) {
    //     CREATEPORE(params.createporego, working_gro, params.poreradius, params.poreaxis)
    //     working_gro = CREATEPORE.out.pore_gro

    //     // update topology
    //     BUILDTOP(params.buildtop_python, working_gro, working_top)
    //     working_top = BUILDTOP.out.output_top
    // }


    // run EM #1
    EM(working_gro, working_top, params.em1)

    // run EM #2 
    EM2(EM.out.em_gro, working_top, params.em2)

    // run EQ #1
    EQ(EM2.out.em_gro, working_top, params.eq1)

    // solvate
    SOLVATE(EQ.out.eq_gro, working_top)

    // create index of System, Solvent, Solute
    MAKEINDEX(SOLVATE.out.solvated_gro)

    // run EM #3
    EM3(SOLVATE.out.solvated_gro, SOLVATE.out.topol, params.em2)

    // run eq #2
    EQ2(EM3.out.em_gro, SOLVATE.out.topol, params.eq1)
 
}
