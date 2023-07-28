#!/usr/bin/env nextflow


params.in = "$baseDir/generate.str"
params.ts2cgLibrary = "$baseDir/Martini3+NLs.LIB"

params.em1 = "$baseDir/mdp/em1.mdp"


/*
 * Run TS2CG on input file
 */
process runTS2CG {
    input:
    path "generate.str"
    path "Martini3+NLs.LIB"

    output:
    file "output.gro"
    
    '''
    PCG -str generate.str  -Bondlength  0.2  -LLIB Martini3+NLs.LIB -function  analytical_shape
    '''
}

process rotate{
    input:
    path "output.gro"

    output:
    file "rotated.gro"

    '''
    gmx editconf -f output.gro -rotate 0 90 0 -o rotated.gro
    '''
}


workflow {
    runTS2CG(params.in, params.ts2cgLibrary) \
      | rotate \
      | view
}