.. image:: docs/source/figures/logo/banner-README.png

.. inclusion-readme-intro-start

NMRforMD (see the `documentation`_ here) is a Python toolkit to calculate NMR relaxation times
from molecular dynamics trajectory files. Used in combination
with `MDAnalysis`_, it allows for the analysis of trajectory
files from any MDAnalysis-compatible simulation package, including
`LAMMPS`_ and `GROMACS`_.

Details about installation, use, and common pitfalls are given in the `documentation`_. 

Notes and known issues
----------------------

- NMRforMD is still in development, please raise an issue here if you encounter a problem
- the code has mostly been tested with GROMACS and LAMMPS trajectory files, but should work with other molecular dynamics packages, as long as they are compatible with MDAnalysis
- NMRforMD does not work with triclinic box, use MDAnalysis to convert your trajectory to orthorhombic
- for very large trajectory file, the code requires a lot of memory
- the code has only beed tested with hydrogen atoms (spin 1/2)
- only works for dipolar interaction, not quadrupolar interaction

.. image:: docs/source/figures/systems/system-white.png

Figure : Example of systems that can be analysed using NMRforMD, from left to right: a 
bulk water reservoir, a PEG molecule, and water in a slit silica pore. 

|mdanalysis| |readthedoc|

.. _`MDAnalysis`: https://www.mdanalysis.org/
.. _`LAMMPS`: https://www.lammps.org/
.. _`GROMACS`: https://www.gromacs.org/
.. _`documentation`: https://nmrformd.readthedocs.io/en/latest/

.. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
    :alt: Powered by MDAnalysis
    :target: https://www.mdanalysis.org

.. |readthedoc| image:: https://readthedocs.org/projects/nmrformd/badge/?version=latest
    :alt: readthedoc
    :target: https://nmrformd.readthedocs.io/en/latest/