/*
 * Build topology
 */
 process BUILDTOP{
    input:
    path "build_top.py"
    path "input.gro"
    path "input.top"

    output:
    path "output.top", emit: output_top

    """
    python build_top.py -c input.gro -p input.top -po output.top
    """
}
