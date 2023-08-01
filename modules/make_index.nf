/*
 * make a solute/solvent index file
 */
process make_index{
    input:
    path "input.gro"

    output:
    path "index.ndx", emit: index

    """
    gmx make_ndx -f input.gro
    """


}