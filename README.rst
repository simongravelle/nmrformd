.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/docs/source/images/NMRforMD_READMEc.png

.. inclusion-readme-intro-start

NMRforMD is a python script to calculate NMR relaxation times from molecular dynamics trajectory files. Used in combination with `MDAnalysis`_, it allows for the analysis of trajectory files from LAMMPS, GROMACS, or AMBER simulation package.

Information
-----------

NMRforMD is in development and likely to return errors. Please raise an issue here if you find one.

.. _`MDAnalysis`: https://www.mdanalysis.org/
.. inclusion-readme-intro-end

For details and a tutorial, have a look at the `documentation`_.
	
Installation
------------

.. inclusion-readme-installation-start

Using pip, type in a terminal:

.. code-block:: bash

	pip install nmrformd

Or, clone this repository on your computer and use pip from the main directory:

.. code-block:: bash

	git clone git@github.com:simongravelle/nmrformd.git
	
	cd nmrformd/

	pip install .
	
You can run the test using pytest:
	
.. code-block:: bash	
	
	cd tests
	pytest mytest.py

.. inclusion-readme-installation-end
.. inclusion-basic-intro-start

Basic example
-------------

This is an example showing how to use NMRforMD to measure NMR signal from 
a molecular dynamics simulations. See the `tutorial`_ for more information.

.. _`tutorial`: https://nmrformd.readthedocs.io/en/latest/documentation_pages/tutorial1.html

.. code-block:: python3

	import MDAnalysis as mda
	import nmrformd
	u = mda.Universe("topology.tpr", "trajectory.xtc")
	nmr_result = nmrformd.NMR(u, "type H", "type H", "full", 0, "m0")

The NMR relaxation time T1 is given by ``nmr_result.T1``.

.. inclusion-basic-intro-end

Known issues
------------

- for very large trajectory file, the code requires a lot of memory
- currently only residues have been tested
- the code has mostly been tested with GROMACS trajectory file

.. _`documentation`: https://nmrformd.readthedocs.io/en/latest/

