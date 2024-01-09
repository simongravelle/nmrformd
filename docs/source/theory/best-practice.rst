Best practice
=============

Choosing the force field
------------------------

.. container:: justify

    The agreement between experiments and simulations can only be as good as the
    force field used in the simulations. Although it has been shown that some
    force fields lead to excellent agreement with experimental data, as for instance
    for water, hydrocarbons, or polymer melts
    :cite:`singerMolecularDynamicsSimulations2017,gravelleNMRInvestigationWater2023,gravelleAssessingValidityNMR2023`,
    it is important to keep in mind that force fields are often parametrized
    to reproduce thermodynamic quantities, such as solvation energy.
    However, NMR relaxation times depend on both structural
    and dynamical quantities, differences between experiments
    and simulations can be expected for the less accurate force field.

.. container:: justify

    Here, as an illustration, the NMR relaxation time :math:`T_1`
    of bulk water was measured as a function of the temperature
    for three different water models:
    :math:`\text{TIP4P}-\epsilon` :cite:`fuentes-azcatlNonPolarizableForceField2014`,
    :math:`\text{SPC/E}` :cite:`berendsenMissingTermEffective1987`,
    and :math:`\text{TIP3P}` :cite:`jorgensenComparisonSimplePotential1983`.
    Our results show that the :math:`\text{TIP4P}-\epsilon` water models
    is in excellent agreement with experimental measurements from 
    Krynicki et al. :cite:`krynickiProtonSpinlatticeRelaxation1966`
    and Hindman et al. :cite:`hindmanRelaxationProcessesWater2003`.
    By contrast, :math:`\text{SPC/E}` and :math:`\text{TIP3P}`
    both overestimate the NMR relaxation time :math:`T_1`, in 
    excellent agreement with previous results
    by Calero et al. :math:`calero1HNuclearSpin2015`. Note that Calero et al.
    used :math:`\text{TIP4P}-2005` water model instead of :math:`\text{TIP4P}-\epsilon`,
    however these two models have very
    similar structures and viscosities :cite:`fuentes-azcatlNonPolarizableForceField2014`.

.. image:: ../figures/illustrations/bulk-water/experimental_comparison-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/illustrations/bulk-water/experimental_comparison-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: figurelegend

    Figure: NMR relaxation time :math:`T_1` between MD simulations of bulk 
    water obtained with three different water models:
    :math:`\text{TIP4P}-\epsilon` :cite:`fuentes-azcatlNonPolarizableForceField2014`,
    :math:`\text{SPC/E}` :cite:`berendsenMissingTermEffective1987`,
    and :math:`\text{TIP3P}` :cite:`jorgensenComparisonSimplePotential1983`.
    Results are compared with experiments 
    from Krynicki et al. :cite:`krynickiProtonSpinlatticeRelaxation1966`
    and from Hindman et al. :cite:`hindmanRelaxationProcessesWater2003`.

Simulation accuracy
-------------------

.. container:: justify

    Since NMR relaxation rate measurements are sensitive both thermodynamic and dynamic quantities, 
    it is important to ensure the accuracy of the simulation.
    For instance, the cut-off for the Lennard-Jones interaction has a slight impact
    on the value :math:`R_1` :cite:`gravelleNMRInvestigationWater2023`.

Box size
--------

.. container:: justify

    NMR relaxation measurements are not extremely sensitive to the box size, however,
    a small effect of the box size can be see, particularly when reaching extremely small boxes:

.. image:: ../figures/best-practices/size-effect-tau-R1-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/size-effect-tau-R1-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: justify

    Looking at the correlation functions, a strong effect of the box size can be 
    see on the inter-molecular contribution, while almost no effect is seen 
    on the intra-molecular contribution:

.. image:: ../figures/best-practices/size-effect-gij-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/size-effect-gij-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: justify

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

.. container:: justify

    For comparison with experimental value, the total duration of the simulation
    must either be larger than :math:`\tau_c`, where :math:`\tau_c` is the longest
    characteristic motion in the system, or be low enough to match the actual Larmor
    frequency used in experiments.

Dumping frequency
-----------------

.. container:: justify

    Dumping period must be smaller than the smaller correlation time of the system, or a 
    significative error on :math:`R_1` will be connected.
