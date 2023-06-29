Bulk water
==========

.. container:: hatnote

   The NMR relaxation time :math:`T_1`  of water is calculated

.. image:: ../../../examples/bulk-water/figures/water-dark-square.png
    :class: only-dark
    :alt: Water molecules simulated with lammps - NMR relaxation time calculation
    :width: 250
    :align: right

.. image:: ../../../examples/bulk-water/figures/water-light-square.png
    :class: only-light
    :alt: Water molecules simulated with lammps - NMR relaxation time calculation
    :width: 250
    :align: right

.. container:: justify

    In this tutorial, the NMR relaxation time :math:`T_1` of water is measured using
    NMRforMD. `MDAnalysis <https://www.mdanalysis.org>`__,
    `numpy <https://www.numpy.org>`__, and
    `matplotlib <https://www.matplotlib.org>`__ and NMRforMD must be
    installed.

    The system is made of 300 TIP4P water molecules simulated in the NPT ensemble with
    `LAMMPS <https://www.lammps.org/>`__ at a temperature of 20°C. The total
    duration of the simulation is 1\,ns, and the timestep is 2 fs. You can
    access the input files in this
    `repository <https://github.com/simongravelle/nmrformd/tree/main/tests>`__,
    which you can use to create larger system or longer trajectory. If
    you are not familiar with LAMMPS, you can find `tutorials
    here <https://lammpstutorials.github.io/>`__.

File preparation
----------------

.. container:: justify

    Either download the files from the Github |repository|, or clone
    the NMRforMD repository:

.. code-block:: bash

    git clone git@github.com:simongravelle/nmrformd.git

.. container:: justify

    The data are located in 'examples/bulk-water/lammps-inputs/'.

    Open a Python script or a Jupyter notebook, and start by defining
    a path to the data files. In my case, since I am working from
    'examples/bulk-water/', its simply 'lammps-inputs/':

.. code-block:: python

	datapath = "lammps-inputs/"

.. |repository| raw:: html

   <a href="ttps://github.com/simongravelle/nmrformd/tree/main/tests" target="_blank">repository</a>

Import the libraries
--------------------

.. container:: justify

    Import numpy, MDAnalysis, and NMRforMD:

.. code-block:: python

	import numpy as np
	import MDAnalysis as mda
	import nmrformd as nmrmd

Create a MDAnalysis universe
----------------------------

.. container:: justify

    Import the configuration file and the trajectory:

.. code-block:: python

	u = mda.Universe(datapath+"topology.data", datapath+"traj.xtc")

.. container:: justify

    The MDAnalysis universe *u* contains both topology (atoms types, masses, etc.)
    and trajectory (atom positions at every frame).

    Let us extract a few information from the universe, such as number of molecules,
    timestep, and total duration:

.. code-block:: python

	n_molecules = u.atoms.n_residues
	print(f"The number of water molecules is {n_molecules}")

>> The number of water molecules is 300

.. code-block:: python

	timestep = np.int32(u.trajectory.dt)
	print(f"The timestep is {timestep} ps")

>> The timestep is 1 ps

.. code-block:: python

	total_time = np.int32(u.trajectory.totaltime)
	print(f"The total simulation time is {total_time} ps")

>> The total simulation time is 1000 ps

Run NMRforMD
------------

..  container:: justify

    Let us isolate a group of atoms containing all the hydrogen atoms (i.e. atoms of 
    type 2) of the system:

.. code-block:: python

	group_i = "type 2"

..  container:: justify

    Then, let us run NMRforMD, using the same group as i and j types:

.. code-block:: python

	nmr_result = nmrmd.NMR(u, [group_i, group_i], number_i=40)

..  container:: justify

    With 'number_i = 40', only 40 randomly selected atoms within 'group_i' are considered for the calculation.
    Increase this number for better resolution. Use 'number_i = 0' to consider all the atoms.

Extract T1
----------

..  container:: justify

    Let us access the calculated value of the NMR relaxation time T1:

.. code-block:: python

	T1 = np.round(nmr_result.T1,2)
	print(f"NMR relaxation time T1 = {T1} s")

>> NMR relaxation time T1 = 3.08 s

The value you get may vary a little, depending on which hydrogen atoms
were randomly selected by NMRforMD.

Plot the spectrum
-----------------

..  container:: justify

    The T1 spectrum can be extracted as 1/nmr_result.R1 (i.e. the invert of R1),
    and the corresponding frequency is given by nmr_result.f. Let up plot
    T1 as a function of f:

.. image:: ../../../examples/bulk-water/figures/T1-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../../../examples/bulk-water/figures/T1-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

Plot the correlation functions
------------------------------

..  container:: justify

    The correlation function Gij can be accessed from nmr_result.gij[0], and the time 
    from nmr_result.t. Let us plot Gij as a function of t:

.. image:: ../../../examples/bulk-water/figures/G-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../../../examples/bulk-water/figures/G-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water