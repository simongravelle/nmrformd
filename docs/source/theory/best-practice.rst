Best practice
=============

Force field
-----------

Force fields are usually parametrized based on thermodynamic quantities only.
However, NMR relaxation quantities depend on both structural and dynamical quantities. 

Simulation precision
--------------------

NMR relaxation rate measurements are extremely sensitive to the precision of the
simulation. The cut-off, for instance, was noted to have a slight impact
on :math:`R_1`.

Box size
--------

NMR relaxation measurements are not extremely sensitive to the box size, however,
a small effect of the box size can be see, particularly when reaching extremely small boxes:

.. image:: ../figures/best-practices/size-effect-tau-R1-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/size-effect-tau-R1-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

Looking at the correlation functions, a strong effect of the box size can be 
see on the inter-molecular contribution, while almost no effect is seen 
on the intra-molecular contribution:

.. image:: ../figures/best-practices/size-effect-gij-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/size-effect-gij-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

Despite the strongly modified correlation functions obtained for small boxes,
the relaxation rate is not so affected:

.. image:: ../figures/best-practices/size-effect-R1-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/size-effect-R1-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

Simulation duration
-------------------

For comparison with experimental value, the total duration of the simulation
must either be larger than :math:`\tau_c`, where :math:`\tau_c` is the longest
characteristic motion in the system, or be low enough to match the actual Larmor
frequency used in experiments.

Dumping frequency
-----------------

Dumping period must be smaller than the smaller correlation time of the system, or a 
significative error on :math:`R_1` will be connected.

For bulk system: check isotropy
-------------------------------

Although generally true for bulk systems, as showed for glycerol :cite:`becherMolecularDynamicsSimulations2021`,
it can be worth ensuring that the relation

.. math::

    \frac{1}{6} G^{(0)} (t) = G^{(1)} (t) = \frac{1}{4} G^{(2)} (t) 

actually stands. For a system of bulk water, the superimposition is clearly verified:

#todo : superimpose water, PEG-water, and slit silica on the same graph

.. image:: ../figures/best-practices/proportionality-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/proportionality-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

If not, all three correlation functions must be calculated.
