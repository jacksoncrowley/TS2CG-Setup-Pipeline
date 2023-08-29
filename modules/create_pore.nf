/*
 * Create Pore
 */
 process CREATEPORE{
    input:
    path "create_pore.go"
    path "input.gro"
    // path "input.top"
    val radius
    val axis

    output:
    path "pore.gro", emit: pore_gro

    script:
    """
    go run create_pore.go -c input.gro -r ${radius} -axis ${axis} -o pore.gro
    """
}
