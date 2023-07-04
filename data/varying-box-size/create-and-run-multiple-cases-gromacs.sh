#!/bin/bash

set -e

gmx=/WORK/simon/gromacs-2023/build-gpu/bin/gmx

for n_mol in 3000 2180 1584 1151 836 608 442 321 233 169 123 89 65 47 34 25 # np.int32(np.logspace(np.log10(25), np.log10(3000), 16))
do
    folder=_N${n_mol}-gromacs/
    if [ ! -d "$folder" ];
    then
        mkdir $folder
    fi
    cd gromacs-inputs/
        newline='nwater = '$n_mol
        oldline=$(cat generate_system.py | grep 'nwater = ')
        sed -i '/'"$oldline"'/c\'"$newline" generate_system.py
    cd ..

    cp -r gromacs-inputs/* $folder

    cd $folder
        ${gmx} grompp -f input/nvt.mdp -o nvt -pp nvt -po nvt
        ${gmx} mdrun -v -deffnm nvt -nt 8 -gpu_id 0 -pin on
        mv nvt.gro conf.gro

        ${gmx} grompp -f input/npt.mdp -o npt -pp npt -po npt
        ${gmx} mdrun -v -deffnm npt -nt 8 -gpu_id 0 -pin on
        mv npt.gro conf.gro

        ${gmx} grompp -f input/run.mdp -o run -pp run -po run
        ${gmx} mdrun -v -deffnm run -nt 8 -gpu_id 0 -pin on
    cd ..
done
