# NMR relaxation time calculator

Python script for the calculation of NMR relaxation time T1 and T2 from molecular dynamics trajectory file. 

### import module

import nmrformd as nmr
import MDAnalysis as mda

### create a universe

u = importuniverse(myfile.tpr,myfile.xtc)

### call nmrtomd

mydic = nmd(u,targetH,neighborH,type,N)

N : number of hydrogen to consired, by default, they will all be considered, otherwise a number N of them will be chosen randomly

type:
- same_res : only consider neighbor hydrogen of the same residue/molecule (for rotational dynamics analysis)
- diff_res : only consider neighbor hydrogen of a different residue/molecule (for translational dynamics analysis)
- all_res : consider all the atoms of neighborH group, except the target hydrogen itself (default)

outputs
- CR
- JR
- CT
- JR
- T1
- T2
- R1 
- T2
- DR
- DT
- DwR
- DwT
