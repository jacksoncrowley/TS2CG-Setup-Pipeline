#!/usr/bin/env nextflow

params.cores = 24


workflow {

    // run TS2CG
    Channel.fromPath(params.in) | runTS2CG

    // rotate TS2CG output
    rotate(runTS2CG.out.output_gro)

    // run em
    // em(rotate.out.rotated_gro, runTS2CG.out.output_top, params.em1)
    view
}


/*
 * Run TS2CG on input file
 */
process runTS2CG {
    publishDir "results/setup", mode:"copy"
    input:
    path input
    // path "Martini3+NLs.LIB"

    output:
    path "output.gro", emit: output_gro
    path "output.top", emit: output_top

    """
    PCG -str ${input} -Bondlength  0.2  -LLIB ${projectDir}/Martini3+NLs.LIB -function  analytical_shape
    rm pcg.log
    """
}

/*
 * rotate a gro file
 */
process rotate{
    publishDir "results/setup", mode:"copy"
    input:
    path "input.gro"

    output:
    path "rotated.gro", emit: rotated_gro

    '''
    gmx editconf -f input.gro -rotate 0 90 0 -o rotated.gro
    '''
}

/*
 * run energy minimisation
 */
process em{
    publishDir "results/setup", mode:"copy"
    input:
    path "input.gro"
    path "input.top"
    path "em.mdp"

    output:
    path "em.gro", emit: em_gro

    """
    gmx grompp -c input.gro -p input.top -f em.mdp -o em.tpr
    gmx mdrun -c em.gro
    """
}
