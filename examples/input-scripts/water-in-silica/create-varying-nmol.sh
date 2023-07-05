#!/bin/bash

set -e

raw_data=../../raw-data/water-in-silica/

# loop over log-spaced data created using np.int32(np.logspace(np.log10(25), np.log10(4000), 12))
for Lz in 50
do

    folder=${raw_data}N${Lz}
    if [ ! -d "$folder" ];
    then
        mkdir $folder
        cd inputs
            newline='Lz = '$Lz
            oldline=$(cat create-system.py | grep 'Lz = ')
            sed -i '/'"$oldline"'/c\'"$newline" create-system.py
            python3 create-system.py
        cd ..
        echo 'Create GROMACS configuration for Lz = '${Lz}
        # Use GROMACS
        cp -r inputs/ff $folder
        cp -r inputs/input ${folder}
        cp inputs/conf.gro $folder
        cp inputs/topol.top $folder
        cp inputs/run_bigfoot_UGA.sh $folder
        rm inputs/conf.gro
        rm inputs/topol.top
    else
        echo 'Folder '${folder}'exist already, skipped'
    fi
done
