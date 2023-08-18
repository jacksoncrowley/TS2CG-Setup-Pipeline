/*
 * rotate a gro file
 */
process ROTATE{
    publishDir "results/setup", mode:"copy"
    input:
    path "input.gro"
    val rotate_str

    output:
    path "rotated.gro", emit: rotated_gro

    """
    gmx editconf -f input.gro -rotate ${rotate_str} -o rotated.gro
    """
}
