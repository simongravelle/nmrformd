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
   
.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/tests/bulk_water_lammps/water-dark.png
    :class: only-dark
    :alt: Water molecules simulated with lammps.

.. image:: https://raw.githubusercontent.com/simongravelle/nmrformd/main/tests/bulk_water_lammps/water-light.png
    :class: only-light
    :alt: Water molecules simulated with lammps.

    Oxygen (type 1) atoms in red, and hydrogen atoms (type 2) in white.
