#!/bin/bash
#OAR -n polymer-water-mixture
#OAR -l /nodes=1/gpu=1/cpu=1/core=8,walltime=48:00:00
#OAR -p gpumodel='A100'
#OAR --stdout emd.out
#OAR --stderr emd.err
#OAR --project tamtam

set -e

export GMX_MAXBACKUP=-1

gmx=/home/gravells/softwares/gromacs-2023/build-gpu/bin/gmx

${gmx} grompp -f input/nvt.mdp -o nvt -pp nvt -po nvt -c conf.gro
${gmx} mdrun -deffnm nvt -v -rdd 1 -nt 8 -pin on

${gmx} grompp -f input/npt.mdp -o npt -pp npt -po npt -c nvt.gro -maxwarn 1
${gmx} mdrun -deffnm npt -v -rdd 1 -nt 8 -pin on

${gmx} grompp -f input/prod.mdp -o prod -pp prod -po prod -c npt.gro -maxwarn 1
${gmx} mdrun -deffnm prod -v -rdd 1 -nt 8 -pin on
