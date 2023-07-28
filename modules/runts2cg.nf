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
