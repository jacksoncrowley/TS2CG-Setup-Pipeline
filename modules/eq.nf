/*
 * run equilibration
 */
process EQ{
    publishDir "${params.outDir}", mode:"copy"
    cpus = params.cores

    input:
    path "input.gro"
    path "input.top"
    path "eq.mdp"

    output:
    path "eq.gro", emit: eq_gro

    """
    gmx grompp -c input.gro -p input.top -f eq.mdp -o eq.tpr -maxwarn 1
    gmx mdrun -s eq.tpr -c eq.gro -nt ${task.cpus} 
    """
}
