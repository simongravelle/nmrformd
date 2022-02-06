Tutorial
========

In this tutorial, the NMR relaxation time T1 of water is going to be measured using NMR for MD. 
`MDAnalysis`_, and NMRforMD must be installed. The data is a short MD trajectory file of bulk
water molecules simulated with `GROMACS`_. If you want to generate longer trajectory files, the 
input files are available in this `repository`_. 

.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/docs/source/images/water_and_T1.png
	:width: 100%

File preparation
----------------

Clone the NMRforMD repository, and go to the tests/bulk_h2o/ folder:

.. code-block:: bash

	git clone git@github.com:simongravelle/nmrformd.git
	cd tests/bulk_h2o/
	
Python script
-------------

Then, in a Python script, import NMRforMD, numpy, and MDAnalysis:

.. code-block:: python3

	import numpy as np
	import nmrformd as NMR
	import MDAnalysis as mda

Create a MDAnalysis universe:

.. code-block:: python3

	u = mda.Universe("topology.tpr","trajectory.xtc")

Let us extract a few information from the universe, just to know what we are dealing with:

.. code-block:: python3

	n_molecules = u.atoms.select_atoms("type OW").atoms.n_atoms
	print(f"The number of water molecules is {n_molecules}")
	timestep = np.round(u.trajectory.dt,2)
	print(f"The timestep is {timestep} ps")
	total_time = np.round(u.trajectory.totaltime,2)
	print(f"The total simulation time is {total_time} ps")

which returns:

.. code-block:: python3

	The number of water molecules is 400
	The timestep is 0.2 ps
	The total simulation time is 100.0 ps

Note that for a proper measurement of T1 for bulk water, a total duration of 2 ns is more reasonable.
Here a smaller trajectory file is used in order to keep the total file size as small as possible.
Let us define two groups for the calculation of the correlation function, see `this paper`_ for details.
Since the only species contributing to the NMR signal in this system are the hydrogen atoms of the water '
molecules, the two groups are identical and both contain all the hydrogen atoms of the system:

.. code-block:: python3

	group_i = "type HW"
	group_j = "type HW"

You can make sure that these groups will be recognized by MDAnalysis:

.. code-block:: python3

	u.atoms.select_atoms(group_i)

which returns:

.. code-block:: python3

	<AtomGroup with 800 atoms>

Then, let us choose a type of analysis, it can be either inter_molecular for the calculation to be 
concerned by the atoms from the same residue (i.e. for rotational dynamics information), intra_molecular
for atoms of different residue (i.e. translational dynamics information), or “full” for both intra and inter.
Let us choose full:

.. code-block:: python3

	analysis = "full"

Choose a number n_i of atom of the group i to consider for the calculation. These atoms will be chosen 
randomly from the group group_i. Here let us use 0 (when 0 is selected, all the atoms of the group are considered which is better here since the statistic is already 
very small):

.. code-block:: python3

	n_i = 0

Finally, choose either "m0" or "m012" for calculation using only the spherical harmonic m=0 (enough for 
isotropic liquid) or all three harmonic m=0, 1 and 2 (can be necessary for more complex system, such 
as water confined in a slit). Here, our system being isotropic, let us choose "m0".

Then, run NMRforMD:

.. code-block:: python3

	nmr_result = NMR.NMR(u, group_i, group_j, analysis, n_i, "m0")

Data analysis
-------------

Let us extract a few quantity from nmr_result, such as the NMR relaxation times T1 and T2 (which are
expected to be equal for water), and the correlation time tau:

.. code-block:: python3
	
	T1 = np.round(nmr_result.T1,2)
	print(f"NMR relaxation time T1 = {T1} s")
	T2 = np.round(nmr_result.T2,2)
	print(f"NMR relaxation time  T2 = {T2} s")
	tau = np.round(nmr_result.tau,2)
	print(f"Correlation time = {tau} ps")

Which returns:

.. code-block:: python3

	NMR relaxation time T1 = 2.58 s
	NMR relaxation time  T2 = 2.58 s
	Correlation time = 3.75 ps

	
The agreement with experiment is not perfect here (experiments give T1 ~ 3s) because of the short 
simulation time. Full relaxation time T1 spectrum can be extracted from 1/nmr_result.R1, and the 
frequency range as nmr_result.f, see the image at the top of this page.

.. _`this paper`: https://www.sciencedirect.com/science/article/abs/pii/S1090780717300319
.. _`MDAnalysis`: https://www.mdanalysis.org
.. _`repository`: https://github.com/simongravelle/nmrformd/tree/main/tests
.. _`GROMACS` : https://www.gromacs.org/

.. toctree::
   :maxdepth: 2
   :caption: NMRforMD
   :hidden:
