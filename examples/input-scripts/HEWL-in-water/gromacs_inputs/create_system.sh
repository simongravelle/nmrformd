#!/bin/bash

export GMX_MAXBACKUP=-1

protein_in_box=protein_in_box.gro
solvated_protein=solvated_protein.gro
final_protein=conf.gro
ion_tpr=ions.tpr

gmx editconf -f reference/protein.gro -o $protein_in_box -c -d 0.5 -bt cubic > /dev/null 2>&1
if test -f "$protein_in_box"
then
	echo "Box size:"
	tail -n 1 $protein_in_box
fi

# change scale to vary the number of water molecules
cp reference/topol.top .
gmx solvate -scale 1.55 -cp $protein_in_box -cs reference/tip4peps.gro -o $solvated_protein -p topol.top > /dev/null 2>&1
# -scale 1.05

echo "Number of molecules:"
tail -n 1 topol.top

#gmx grompp -f input/ions.mdp -c $solvated_protein -p topol.top -o $ion_tpr > /dev/null 2>&1
#gmx genion -s $ion_tpr -o $final_protein -p topol.top -pname NA -nname CL -neutral <<EOF > /dev/null 2>&1
#    13
#EOF

gmx insert-molecules -f $solvated_protein -ci reference/CL.gro -o conf.gro -nmol 8 -try 500 > /dev/null 2>&1
#cp $solvated_protein conf.gro
echo "CL 8" >> topol.top

echo "Number of ions:"
tail -n 1 topol.top

for T in 300
do
	folder="../../../raw-data/HEWL-in-water/T"${T}"K/"

	if [ ! -d "$folder" ];
  	then
  		echo "Creating folder T"${T}"K"
  		mkdir ${folder}
  	fi

	newline='ref-t = '${T}' '${T}' '${T}' ; reference temperature for coupling (K)'
    oldline=$(cat input/run.mdp  | grep 'ref-t = ')
    sed -i '/'"$oldline"'/c\'"$newline" input/run.mdp 
    
    newline='ref-t = '${T}' '${T}' '${T}' ; reference temperature for coupling (K)'
    oldline=$(cat input/nvt.mdp  | grep 'ref-t = ')
    sed -i '/'"$oldline"'/c\'"$newline" input/nvt.mdp
  
    newline='gen-temp = '${T}' ; temperature for Maxwell distribution (K)'
    oldline=$(cat input/nvt.mdp  | grep 'gen-temp = ')
    sed -i '/'"$oldline"'/c\'"$newline" input/nvt.mdp
   
  	cp -r perso.oplsaa.ff/ ${folder}perso.oplsaa.ff
  	cp -r input/ ${folder}input
  	cp run_bigfoot_UGA.sh ${folder}run_bigfoot_UGA.sh
  	cp conf.gro ${folder}conf.gro
  	cp topol.top ${folder}topol.top
  	cp evaluate_mass_ratio.py ${folder}evaluate_mass_ratio.py

    python3 evaluate_mass_ratio.py
done


















