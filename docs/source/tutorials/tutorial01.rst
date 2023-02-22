Tutorial 1 : bulk water
=======================

Here the NMR relaxation time :math:`T_1` of water is measured using
NMRforMD. `MDAnalysis <https://www.mdanalysis.org>`__,
`numpy <https://www.numpy.org>`__, and
`matplotlib <https://www.matplotlib.org>`__ and NMRforMD must be
installed.

The system is a short molecular dynamics trajectory of bulk TIP4P
water molecules simulated in the NPT ensemble simulated with
`LAMMPS <https://www.lammps.org/>`__ (temperature 20°C). You can
access the input files in this
`repository <https://github.com/simongravelle/nmrformd/tree/main/tests>`__,
which you can use to create larger system or longer trajectory. If
you are not familiar with LAMMPS, you can find `tutorials
here <https://lammpstutorials.github.io/>`__.
   
.. image:: ../../../tests/bulk_water_lammps/water-dark.png
    :class: only-dark
    :alt: Water molecules simulated with lammps.

.. image:: ../../../tests/bulk_water_lammps/water-light.png
    :class: only-light
    :alt: Water molecules simulated with lammps.

Oxygen (type 1) atoms in red, and hydrogen atoms (type 2) in white.

File preparation
----------------

Either download the files from the Github [repository](https://github.com/simongravelle/nmrformd/tree/main/tests), or clone the NMRforMD repository:

.. code-block:: bash

    git clone git@github.com:simongravelle/nmrformd.git

The data are located in 'tests/bulk_water_lammps/'.

Open a Python script or a Jupyter notebook, and start by defining a path to the data files. In
my case, its:

.. code-block:: python

	datapath = "tests/bulk_water_lammps/"

Import the libraries
--------------------

Import numpy, MDAnalysis, pyplot, and NMRforMD:

.. code-block:: python

	import numpy as np
	import MDAnalysis as mda
	import nmrformd as nmrmd

Create a MDAnalysis universe
----------------------------

Import the configuration file and the trajectory:

.. code-block:: python

	u = mda.Universe(datapath+"conf.data", datapath+"traj.xtc")

The MDAnalysis universe *u* contains both topology (atoms types, masses, etc.)
and trajectory (atom positions every frame).

Let us extract a few information from the universe:

.. code-block:: python

	n_molecules = u.atoms.select_atoms("type 1").atoms.n_atoms
	print(f"The number of water molecules is {n_molecules}")

.. code-block:: bash

	>> The number of water molecules is 400

.. code-block:: python

	timestep = np.int32(u.trajectory.dt)
	print(f"The timestep is {timestep} ps")

.. code-block:: bash

	>> The timestep is 1 ps

.. code-block:: python

	total_time = np.int32(u.trajectory.totaltime)
	print(f"The total simulation time is {total_time} ps")

.. code-block:: bash

	>> The total simulation time is 1000 ps

Run NMRforMD
------------

Let us isolate a group of atoms containing all the hydrogen atoms of the system:

.. code-block:: python

	group_i = "type 2"

Then, let us run NMRforMD:

.. code-block:: python

	nmr_result = nmrmd.NMR(u, [group_i, group_i], number_i=40)

With 'number_i = 40', only 40 randomly selected atoms within 'group_i' are considered for the calculation.
Increase this number for better resolution. Use 'number_i = 0' to consider all the atoms.

Extract T1/T2
-------------

Access the calculated value of NMR relaxation time T1:

.. code-block:: python

	T1 = np.round(nmr_result.T1,2)
	print(f"NMR relaxation time T1 = {T1} s")

.. code-block:: bash

	>> NMR relaxation time T1 = 3.08 s

Plot the spectrum and the correlation functions
-----------------------------------------------

.. code-block:: python

	from matplotlib import pyplot as plt

	fig = plt.figure(figsize=(14,7))
	ax1 = plt.subplot(1, 2, 1)
	ax1.loglog(nmr_result.f[:-250], 1/nmr_result.R1[:-250], 'o', markersize=8)
	ax1.set_xlabel(r"$f$ (MHz)", fontdict=font)
	ax1.set_ylabel(r'$T_1$ (s)', fontdict=font)
	plt.xlim(5e2, 5e5)
	plt.ylim(1, 100)
	ax2 = plt.subplot(1, 2, 2)
	ax2.semilogx(nmr_result.t[:-250], nmr_result.gij[0][:-250], 'o', markersize=8)
	ax2.set_xlabel(r"$t$ (ps)", fontdict=font)
	ax2.set_ylabel(r'$C$', fontdict=font)
	plt.xlim(5e-1, 5e2)
	plt.ylim(-0.5e10, 5e10)
	plt.show()

.. image:: ../../../tests/bulk_water_lammps/spectrums-dark.png
    :class: only-dark
    :alt: NMR results obtained from the simulation.

.. image:: ../../../tests/bulk_water_lammps/spectrums-light.png
    :class: only-light
    :alt: NMR results obtained from the simulation.