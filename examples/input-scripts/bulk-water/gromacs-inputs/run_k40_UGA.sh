#!/bin/bash

set -e

gmx=/WORK/simon/gromacs-2023/build-gpu/bin/gmx

${gmx} grompp -f input/nvt.mdp -c conf.gro -o nvt -pp nvt -po nvt
${gmx} mdrun -v -deffnm nvt -nt 8 -gpu_id 0 -pin on

${gmx} grompp -f input/npt.mdp -c nvt.gro -o npt -pp npt -po npt
${gmx} mdrun -v -deffnm npt -nt 8 -gpu_id 0 -pin on

${gmx} grompp -f input/run.mdp -c npt.gro -o run -pp run -po run
${gmx} mdrun -v -deffnm run -nt 8 -gpu_id 0 -pin on

