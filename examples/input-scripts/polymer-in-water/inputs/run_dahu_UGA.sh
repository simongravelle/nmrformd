#!/bin/bash
#OAR -n NMR-test
#OAR -l /nodes=1/cpu=1/core=8,walltime=48:00:00
#OAR --stdout log.out
#OAR --stderr log.err
#OAR --project tamtam

lmp=/home/gravells/softwares/lammps-8Feb2023/src/lmp_mpi

mpirun -np 8 ${lmp} -in input.lammps
