/*
 * Run TS2CG on input file
 */
process RUNTS2CG {
    input:
    path pcg
    path input
    path top_header

    output:
    path "output.gro", emit: output_gro
    path "system.top", emit: system_top

    """
    ./${pcg} -str ${input} -Bondlength  0.2  -LLIB ${baseDir}/top/Martini3+NLs.LIB -function  analytical_shape
    rm pcg.log

    cat ${top_header} > system.top
    cat output.top >> system.top

    echo ${baseDir}
    sed -i 's:BASEDIR:${baseDir}:g' system.top
    """
}
