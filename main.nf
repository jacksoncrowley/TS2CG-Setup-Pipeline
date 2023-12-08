#!/usr/bin/env nextflow

// paths to required input files
params.pcg          = "$baseDir/PCG"            // PCG from TS2CG
params.em1          = "$baseDir/mdp/em1.mdp"
params.em2          = "$baseDir/mdp/em2.mdp"
params.eq1          = "$baseDir/mdp/eq1.mdp"
params.top_header   = "$baseDir/top/header.txt" // file listing .itps in order
params.input        = "$launchDir/generate.str"   // input .str file
params.outDir       = "$launchDir/results"

// number of processes to be used by energy minimization and equilibration
params.cores = 16

// import modules
include { RUNTS2CG } from './modules/runts2cg.nf'
include { EM; EM as EM2; EM as EM3 } from './modules/em.nf'
include { EQ; EQ as EQ2 } from './modules/eq.nf'
include { SOLVATE } from './modules/solvate.nf'
include { MAKEINDEX } from './modules/make_index.nf'


workflow {
    // run TS2CG
    RUNTS2CG(params.pcg, params.input, params.top_header)
    working_gro = RUNTS2CG.out.output_gro
    working_top = RUNTS2CG.out.system_top

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
