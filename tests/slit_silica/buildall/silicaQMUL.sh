#!/bin/sh
#$ -cwd
#$ -pe smp 8        # 8 cores per single GPU (32 per gpu node)
#$ -l h_vmem=7.5G   # We recommend requesting 7.5G per core
#$ -l gpu=1
#$ -l h_rt=240:0:0   # Generally we suggest to use 240 hours

module unload cuda
module load gromacs/5.1.4-gpu

gmx grompp -f ../input/nvt.mdp -o nvt -pp nvt -po nvt -r conf.gro
gmx mdrun -ntomp ${NSLOTS} -v -deffnm nvt -rdd 0.8
mv nvt.gro conf.gro

gmx grompp -f ../input/npt1000.mdp -o npt1000 -pp npt1000 -po npt1000 -r conf.gro
gmx mdrun -ntomp ${NSLOTS} -v -deffnm npt1000 -rdd 0.8
mv npt1000.gro conf.gro

gmx grompp -f ../input/npt.mdp -o npt -pp npt -po npt -r conf.gro
gmx mdrun -ntomp ${NSLOTS} -v -deffnm npt
mv npt.gro conf.gro

gmx grompp -f ../input/runHR.mdp -o runHR -pp runHR -po runHR -r conf.gro
gmx mdrun -ntomp ${NSLOTS} -v -deffnm runHR
mv runHR.gro conf.gro

#gmx grompp -f ../input/runLR.mdp -o runLR -pp runLR -po runLR -r conf.gro
#gmx mdrun -ntomp ${NSLOTS} -v -deffnm runLR

#rm *#*

