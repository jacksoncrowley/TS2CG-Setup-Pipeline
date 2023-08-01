/*
 * solvate (and neutralise?) AND generate new index files
 */
process solvate{
    publishDir "results/setup", mode:"copy"
    input:
    path "input.gro"
    path "topol.top"
    path "water.gro"

    output:
    path "solvated.gro", emit: solvated_gro
    path "topol.top", emit: topol

    """
    gmx solvate -cs water.gro -cp input.gro -p topol.top -radius 0.21 -o solvated.gro
    """
}
