#!/bin/bash

# for smaller boxes, LAMMPS is used due to some GROMACS restriction
# for larger boxes, GROMACS is used as its faster

set -e

raw_data=../../raw-data/bulk-water/

# loop over log-spaced data created using np.int32(np.logspace(np.log10(25), np.log10(4000), 12))
for n_mol in 25 39 62 99 158 251 398 631 1002 1589 2521 4000
do

    folder=${raw_data}N${n_mol}
    if [ ! -d "$folder" ];
    then
        mkdir $folder
        if ((${n_mol} > 700));
        then
            echo 'Create GROMACS configuration for Nmol = '${n_mol}
            # Use GROMACS
            cd gromacs-inputs/
                newline='nwater = '$n_mol
                oldline=$(cat generate_system.py | grep 'nwater = ')
                sed -i '/'"$oldline"'/c\'"$newline" generate_system.py
    	        python3 generate_system.py
             cd .. 
             cp -r gromacs-inputs/ff $folder
             cp -r gromacs-inputs/input $folder
             cp gromacs-inputs/conf.gro $folder
             cp gromacs-inputs/topol.top $folder
             cp gromacs-inputs/run_k40_UGA.sh $folder
             rm gromacs-inputs/conf.gro
             rm gromacs-inputs/topol.top
        else
            echo 'Create LAMMPS configuration for Nmol = '${n_mol}
            # Use LAMMPS
            cd lammps-inputs/
                newline='variable nml equal '$n_mol
                oldline=$(cat input.lammps | grep 'variable nml equal')
                sed -i '/'"$oldline"'/c\'"$newline" input.lammps
            cd ..
            cp lammps-inputs/*.lammps $folder
            cp lammps-inputs/*.sh $folder 
        fi
    else
        echo 'Folder '${folder}'exist already, skipped'
    fi
done
