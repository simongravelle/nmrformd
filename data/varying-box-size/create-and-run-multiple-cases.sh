#!/bin/bash

set -e

lmp=/home/simon/Softwares/lammps-28Mar2023/src/lmp_mpi

for n_mol in 49 89 161 289 518 931 1671 3000 #  np.int32(np.logspace(np.log10(50), np.log10(3000), 8))
do
    folder=_N${n_mol}/
    if [ ! -d "$folder" ];
    then
        mkdir $folder
    fi
    cp lammps-inputs/* $folder

    cd $folder
    mpirun -np 8 ${lmp} -in input.create.lammps -var nml ${n_mol}
    mpirun -np 8 ${lmp} -in input.prod.lammps
    cd ..
done