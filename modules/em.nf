/*
 * run energy minimisation
 */
process em{
    publishDir "results/setup", mode:"copy"
    input:
    path "input.gro"
    path "input.top"
    path "em.mdp"

    output:
    path "em.gro", emit: em_gro

    """
    gmx grompp -c input.gro -p input.top -f em.mdp -o em.tpr -maxwarn 1
    gmx mdrun -s em.tpr -c em.gro -nt 24 
    """
}
