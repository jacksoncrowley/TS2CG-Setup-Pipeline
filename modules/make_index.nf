/*
 * make a solute/solvent index file
 */
process MAKEINDEX{
    publishDir "${params.outDir}", mode:"copy"
    input:
    path "input.gro"

    output:
    path "index.ndx", emit: index

    """
    echo -e "del 1-99\na W NA CL CA\n0 & ! 1\nname 1 Solvent\nname 2 Solute\nq" | gmx make_ndx -f input.gro
    """
}