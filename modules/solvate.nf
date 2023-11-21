/*
 * solvate (and neutralise?) AND generate new index files
 */
process SOLVATE{
    publishDir "${params.outDir}", mode:"copy"
    input:
    path "input.gro"
    path "topol.top"

    output:
    path "solvated.gro"
    path "topol.top", emit: topol

    """
    gmx solvate -cs ${baseDir}/top/water.gro -cp input.gro -p topol.top -radius 0.21 -o solvated.gro
    """
}
