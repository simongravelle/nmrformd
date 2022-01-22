NMRforMD
########

NMRforMD is a python script to calculate NMR relaxation times from molecular dynamics trajectory files. Used in combination with MDAnalysis, it allows for the analysis of trajectory files from LAMMPS, GROMACS, or AMBER simulation package.

Information
-----------

NMRforMD is in development and likely to return errors. Please raise an issue here if you find one. 

How to
------

In a python script, import NMRforMD and MDAnalysis:

.. code-block:: python3

	import nmrformd as NMR
	import MDAnalysis as mda

Then create a MDAnalysis universe:

.. code-block:: python3

	u = importuniverse(topology.tpr,trajectory.xtc)

You can use the topology and trajectory files given in the tests folder. Then, choose the type of atoms for the correlation. Two groups are needed, they can be different, or the same, depending of what you want to achieve: 

.. code-block:: python3

	group_i = "type HW"
	group_j = group_i

Choose a type of analysis, it can be either "inter_molecular" if you only want the calculation to be concerned by the atoms from the same residue (i.e. for rotational dynamics information), "intra_molecular" for atoms of different residue (i.e. translational dynamics information), or "full" for both intra and inter.

.. code-block:: python3

	analysis = "full"

Choose a number of atom of the group i to consider for the calculation. 

.. code-block:: python3
	
	n_i = 100

Finally, choose either "m0" or "m012" for calculation using only the spherical harmonic m=0 (enough for isotropic liquid) or all three harmonic m=0, 1 and 2 (can be necessary for more complex system). Then, run NMRforMD:

.. code-block:: python3

	nmr_result = NMR.NMR(u, groupe_i, group_j, analysis, n_i, "m0")

Results can be accessed from the nmr_result. 
