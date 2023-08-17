/*
 * solvate (and neutralise?) AND generate new index files
 */
process solvate{
    publishDir "results/setup", mode:"copy"
    input:
    path "input.gro"
    path "topol.top"

    output:
    path "solvated.gro", emit: solvated_gro
    path "topol.top", emit: topol

    """
    gmx solvate -cs ${baseDir}/top/water.gro -cp input.gro -p topol.top -radius 0.21 -o solvated.gro
    """
}
