.. image:: https://github.com/simongravelle/nmrformd/blob/main/docs/source/images/NMRforMD_READMEc.png

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

.. inclusion-readme-installation-end

Known issues
------------

- for very large trajectory file, the code requires a lot of memory
- currently only residue are accepted to differentiate atoms from the same molecule/structure/residue
- the code has mostly be tested with GROMACS trajectory file

.. _`documentation`: https://nmrformd.readthedocs.io/en/latest/

