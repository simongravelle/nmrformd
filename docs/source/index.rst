NMR fom MD
==========

.. container:: justify

    Dipolar Nuclear Magnetic Resonance for Molecular Dynamics
    (NMRforMD or formerly NMRforMD) simulations
    is a Python toolkit designed for the computation of
    dipolar NMR relaxation times
    (the so called :math:`T_1` and :math:`T_2`)
    from molecular dynamics simulations.
    Used in combination with |MDAnalysis|,
    NMRforMD allows for the analysis of trajectory
    files from any MDAnalysis-compatible simulation package, including
    |LAMMPS| and |GROMACS|.

.. |MDAnalysis| raw:: html

   <a href="https://www.mdanalysis.org/" target="_blank">MDAnalysis</a>

.. |LAMMPS| raw:: html

   <a href="https://www.lammps.org/" target="_blank">LAMMPS</a>

.. |GROMACS| raw:: html

   <a href="https://www.gromacs.org/" target="_blank">GROMACS</a>

.. image:: ../../avatars/avatars.png
    :class: only-dark
    :alt: molecular dynamics systems used in these examples 

.. image:: ../../avatars/avatars.png
    :class: only-light
    :alt: molecular dynamics systems used in these examples

.. container:: figurelegend

    Figure: Examples of systems that can be analysed
    using NMRforMD, from left to right: a bulk water system, 
    a lennard-jones fluid, and a lysozyme in water.

Datasets
--------

.. container:: justify

    Two molecular dynamics datasets are available on Github. One 
    is a |polymer in water| system generated using LAMMPS, 
    the second is a |water confined in silica|
    generated using GROMACS. The datasets are
    provided to follow the tutorials.

.. |polymer in water| raw:: html

   <a href="https://github.com/simongravelle/polymer-in-water.git" target="_blank">polymer in water</a>

.. |water confined in silica| raw:: html

   <a href="https://github.com/simongravelle/water-in-silica.git" target="_blank">water confined in silica</a>

.. toctree::
   :maxdepth: 2
   :hidden:

   modules/NMR
   modules/utilities

.. toctree::
   :maxdepth: 2
   :caption: Tutorials
   :hidden:
   
   tutorials/installation
   tutorials/isotropic-system
   tutorials/anisotropic-system

.. toctree::
   :maxdepth: 2
   :caption: Illustrations
   :hidden:

   illustrations/lennard-jones-fluids
   illustrations/bulk-water
   illustrations/lysozyme-in-water

.. toctree::
   :maxdepth: 2
   :caption: Theory
   :hidden:

   theory/context
   theory/dipolar_relaxation
   theory/best-practice

.. toctree::
   :maxdepth: 2
   :caption: Additional
   :hidden:

   additional/bibliography
   additional/acknowledgments
