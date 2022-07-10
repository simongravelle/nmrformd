.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/docs/source/images/NMRforMD_READMEc.png

.. inclusion-readme-intro-start

NMRforMD is a python toolkit to calculate NMR relaxation times
from molecular dynamics trajectory files. Used in combination
with `MDAnalysis`_, it allows for the analysis of trajectory
files from `LAMMPS`_ and `GROMACS`_ simulation package.

Information
-----------

This documentation is separated in four parts: tutorials, how-to scripts,
description, and theory.

.. _`MDAnalysis`: https://www.mdanalysis.org/
.. _`LAMMPS`: https://www.lammps.org/
.. _`GROMACS`: https://www.gromacs.org/
.. inclusion-readme-intro-end

For details and instructions for beginners,
have a look at the `documentation`_.

Notes :
    - NMRforMD is still in development, please raise an issue here if you encounter a problem.
    - the code has mostly been tested with GROMACS and LAMMPS trajectory files, but should work with other molecular dynamics packages

Installation
------------

.. inclusion-readme-installation-start

Using pip, type in a terminal:

.. code-block:: bash

	pip3 install nmrformd

.. inclusion-readme-installation-end

To get the last version, clone this repository on your computer
and use pip3 from the main directory:

.. code-block:: bash

	git clone https://github.com/simongravelle/nmrformd.git
	
	cd nmrformd/

	pip install .
	
You can run the tests using pytest:
	
.. code-block:: bash	
	
	cd tests
	pytest .

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
	nmr_result = nmrformd.NMR(u, "type H", "type H")

The NMR relaxation time T1 is given by ``nmr_result.T1``.

.. inclusion-basic-intro-end

Known issues
------------

- for very large trajectory file, the code requires a lot of memory

.. _`documentation`: https://nmrformd.readthedocs.io/en/latest/

