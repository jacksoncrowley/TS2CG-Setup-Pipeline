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
