#!/bin/bash

set -e

jupyter nbconvert --to script 'build_system.ipynb'

for surface in 4;
do
	newline='id_silica = '$surface
	linetoreplace=$(cat build_system.py | grep 'id_silica =')
	sed -i '/'"$linetoreplace"'/c\'"$newline" build_system.py

 	DIRNAME=Id_silica_${surface}
 	#DIRNAME=N4000
	if [ ! -d "../$DIRNAME" ]; then
		echo "creating folder"
		mkdir ../$DIRNAME	
	fi

	python3 build_system.py

	cp posre_SiOH.itp silica.itp ff/
	cp -r ff conf.gro topol.top silicaQMUL.sh ../$DIRNAME
done


