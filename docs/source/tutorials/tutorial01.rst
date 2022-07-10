Tutorial
========

In this tutorial, the NMR relaxation time T1 of water is going to be measured using 'NMRforMD'. 
`MDAnalysis`_, `numpy`_, and `matplotlib`_ and NMRforMD must be installed. 

The system is a short molecular dynamics trajectory of bulk TIP4P water molecules in the NVT ensemble 
simulated with `LAMMPS`_ (temperature 20°C). 

If you want to generate longer trajectory files, the 
input files are available in this `repository`_. 

.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/docs/source/images/water_and_T1.png
	:width: 100%

File preparation
----------------

Clone the NMRforMD repository, and go to the tests/bulk_water_lammps_large/ folder:

.. code-block:: bash

	git clone git@github.com:simongravelle/nmrformd.git
	cd tests/bulk_h2o/
	
Import the libraries
--------------------

Open a Python script, import numpy, MDAnalysis, pyplot, and NMRforMD:

.. code-block:: python3

	import numpy as np
	import MDAnalysis as mda
	from matplotlib import pyplot as plt
	import nmrformd as nmrmd

Create a MDAnalysis universe
----------------------------

Import the configuration file and the trajectory:

.. code-block:: python3

	u = mda.Universe("conf.data", "traj.xtc")

Let us extract a few information from the universe:

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
	The timestep is 1.0 ps
	The total simulation time is 1000.0 ps

Run NMRforMD
------------

Let us isolate a group of atoms containing all the hydrogen atoms of the system:

.. code-block:: python3

	group_i = "type 2"

Then, let us run NMRforMD:

.. code-block:: python3

	nmr_result = nmrmd.NMR(u, [group_i, group_i], number_i = 40)

With 'number_i = 40', only 40 atoms within 'group_i' are considered for the calculation. 
Increase this number for better resolution. Use 'number_i = 0' to consider all the atoms.

Extract T1/T2
-------------

Access the calculated value of T1 and T2 using:

.. code-block:: python3
	
	T1 = np.round(nmr_result.T1,2)
	print(f"NMR relaxation time T1 = {T1} s")
	T2 = np.round(nmr_result.T2,2)
	print(f"NMR relaxation time  T2 = {T2} s")

which returns:

.. code-block:: python3

	NMR relaxation time T1 = 2.81 s
	NMR relaxation time T2 = 2.81 s

Note that for a bulk water system, T1 is known to be equal to T2.

Plot the spectrum and the correlation functions
-----------------------------------------------

.. code-block:: python3

	fig = plt.figure(figsize=(10,5))
	ax1 = plt.subplot(1, 2, 1)
	ax1.loglog(nmr_result.f[:-250], nmr_result.R1[:-250], '.')
	ax1.set_xlabel(r"$f$ (MHz)")
	ax1.set_ylabel(r'$R_1$ (s$^{-1}$)')
	ax2 = plt.subplot(1, 2, 2)
	ax2.semilogx(nmr_result.t[:-500], nmr_result.gij[0][:-500], '.')
	ax2.set_xlabel(r"$t$ (ps)")
	ax2.set_ylabel(r'$C$')
	fig.tight_layout()
	plt.show()

.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/docs/source/images/bulk_water_R1_C.png
	:width: 100%

.. _`this paper`: https://www.sciencedirect.com/science/article/abs/pii/S1090780717300319
.. _`MDAnalysis`: https://www.mdanalysis.org
.. _`numpy`: https://www.numpy.org
.. _`matplotlib`: https://www.matplotlib.org
.. _`repository`: https://github.com/simongravelle/nmrformd/tree/main/tests
.. _`LAMMPS`: https://www.lammps.org/

.. toctree::
   :maxdepth: 2
   :caption: NMRforMD
   :hidden:
