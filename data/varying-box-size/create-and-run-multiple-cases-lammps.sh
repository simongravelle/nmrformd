#!/bin/bash

set -e

lmp=/home/simon/Softwares/lammps-28Mar2023/src/lmp_mpi

for n_mol in 25 34 47 65 89 123 169 233 321 442 608 836 1151 1584 2180 3000 # np.int32(np.logspace(np.log10(25), np.log10(3000), 16))
do
    folder=_N${n_mol}-lammps/
    if [ ! -d "$folder" ];
    then
        mkdir $folder
    fi
    cp lammps-inputs/* $folder

    cd $folder
        newline='variable nml equal '$n_mol
        oldline=$(cat input.create.lammps | grep 'variable nml equal')
        sed -i '/'"$oldline"'/c\'"$newline" input.create.lammps
        mpirun -np 8 ${lmp} -in input.create.lammps
        mpirun -np 8 ${lmp} -in input.prod.lammps
    cd ..
done
