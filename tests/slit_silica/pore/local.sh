#!/bin/bash

gmx grompp -f ../input/prod.mdp -o prod -pp prod -po prod -r conf.gro
gmx mdrun -v -deffnm prod
mv prod.gro conf.gro
