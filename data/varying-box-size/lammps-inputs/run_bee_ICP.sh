#!/bin/bash
#!/bin/sh 
#SBATCH --job-name=NMR-test
#SBATCH --mail-user=sgravelle@icp.uni-stuttgart.de
#SBATCH --mail-type=FAIL
#SBATCH -n 8
#SBATCH --time=48:00:0
#SBATCH --mem=6000

lmp=/home/sgravelle/lammps-23Jun2022/src/lmp_mpi

module load gcc/10.2.0
module load openmpi

srun -n 8 ${lmp} -in input.create.lammps
srun -n 8 ${lmp} -in input.prod.lammps