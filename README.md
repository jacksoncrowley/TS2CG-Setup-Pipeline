# TS2CG Setup Pipeline
This is a basic pipeline for the program [TS2CG](https://github.com/marrink-lab/TS2CG1.1). Further reading can be found at (https://doi.org/10.1038/s41467-020-16094-y). Please cite this paper in case of any starting structures which are generated from this workflow for eventual publication use.

More specifically, the program takes a .str file for the membrane builder program, PCG, and performs the required energy minimisation, solvation and equilibration stages. For further information, see the [Martini 2021 Online Workshop Tutorial](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/558-9-ts2cg), starting at "Tut6: Fixed shapes".

**This workflow currently requires that you have the PCG program installed in your $PATH directory**

### Prerequisites
- Gromacs
- Nextflow
- [TS2CG PCG program](https://github.com/marrink-lab/TS2CG1.1/blob/master/PCG) (included by default in base directory)

## Usage
Within "main.nf", the following variables are defined:

~~~
params.in = "$baseDir/generate.str"
params.outDir = "$baseDir/results"
params.em1 = "$baseDir/mdp/em1.mdp"
params.em2 = "$baseDir/mdp/em2.mdp"
params.eq1 = "$baseDir/mdp/eq1.mdp"
params.top_header = "$baseDir/top/header.txt"
~~~

At a minimum, the user needs only to change the content of, and perhaps the name of generate.str; the PCG input file for generation of the structure.

Providing an writeout directory via the command line can be done by the --outDir argument. In Nextflow, the absolute path must be provided.

Using the argument --rotate, users can specify an optional rotation of the system by providing the argument as degrees in the format

~~~
--rotate x,y,z
~~~
in practice;
~~~
--rotate 0,90,0
~~~
To rotate the system by 90 degrees in y.

Next most likely would be the file "top/header.txt", which is simply the #include statements for the gromacs .itp files to be added to the top of the .top file produced by TS2CG. Indeed, there is a smarter way to do this; expect in a future minor update.

### GROMACS version notice
The first energy minimisation conducted, mdp/em1.mdp, uses soft core potentials, which can produce a bug in certain versions of GROMACS. 

Using GROMACS/2022.5 yielded the following error:
>   Fatal error:
  There are perturbed non-bonded pair interactions beyond the pair-list cutoff
  of 1.155 nm, which is not supported. This can happen because the system is
  unstable or because intra-molecular interactions at long distances are
  excluded. If the latter is the case, you can try to increase nstlist or rlist
  to avoid this.The error is likely triggered by the use of couple-intramol=no
  and the maximal distance in the decoupled molecule exceeding rlist.

The same error was not found using GROMACS/2019.6, GROMACS/2020.7, or GROMACS/2023.1.