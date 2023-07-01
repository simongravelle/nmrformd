
Theory
======

Let us consider an ensemble of identical spins, characterized by the gyromagnetic
ratio :math:`\gamma` and the spin quantum number :math:`I`. The magnetic dipolar
interaction between two spins, :math:`i` and :math:`j`, may be written in terms
of the Hamiltonian :cite:`grivetNMRRelaxationParameters2005` :cite:`bloembergenRelaxationEffectsNuclear1948`:

.. math::

    \hat H_d = \dfrac{\mu_0}{4 \pi} \hbar^2 \gamma^2 \sum_{-2}^{2} F_{ij}^{(m)} (t) \hat A_{ij}^{(m)},

where :math:`\hat A_{ij}^{(m)}` are dimensionless spin operators involving spins :math:`i` and :math:`j`,
and :math:`F_{ij}^{(m)} (t)` are functions of the vector :math:`\boldsymbol{r}_{ij}` between
spins :math:`i` and :math:`j`, which varies through time due to molecular motion.
The function :math:`{\cal F}_{ij} (t)` reads

.. math::
    
    {\cal F}_{ij} (t) = \alpha_m \dfrac{1}{r_{ij}^3 (t)} Y^{(m)}_2 (\Omega_{ij}),

where :math:`Y^{(m)}_2` are normalised spherical harmonics, and where


.. math::

    \alpha_0 = \sqrt{\frac{16 \pi}{5}}, ~ \alpha_1 = \sqrt{\frac{8 \pi}{15}}, ~ \alpha_2 = \sqrt{\frac{32 \pi}{15}}.

:math:`\Omega_{ij}` denotes the polar angles of the direction of :math:`\boldsymbol{r}_{ij}` with respect
to laboratory axes, assuming that the applied static magnetic field is parallel to :math:`\boldsymbol{e}_z`.

Relaxation rate calculation relies on the evaluation of the correlation functions

.. math::

    G^{(m)} (t) = \dfrac{1}{N}
    \sum_{i \ne j}^{N} \left< {\cal F}_{ij}^{(m)} (0) {\cal F}_{ij}^{(m)} (t)  \right>,

where :math:`N` is the number of spin pairs. Spectral densities are obtained from the
Fourier transforms of the correlation functions, 

.. math::

    J^{(m)} (\omega) = \int_\infty^\infty G^{(m)} (t) \mathrm e^{- i \omega t} \mathrm dt 

from which the relaxation rates can be calculated as

.. math::

    R_1 &=&  K \left[ J^{(1)} (\omega_0) + J^{(2)} (2 \omega_0) \right],

    R_2 &=& \dfrac{1}{4} K \left[ J^{(0)} (0) + 10 J^{(1)} (\omega_0) + J^{(2)} (2 \omega_0) \right],

were :math:`\omega_0 = \gamma B_0` is the Larmor frequency, :math:`I = 1/2` the
spin quantum number, and

.. math::

    K = \dfrac{3}{2}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

where :math:`\mu_0` is the vacuum permeability.

Separating intra-molecular and inter-molecular
----------------------------------------------

Intra-molecular and inter-molecular contributions to :math:`R_1`
can be extracted separately, by splitting the correlation functions as:

.. math::

    G^{(m)}_\text{R, T} (t) = \dfrac{1}{N_\text{R, T}}
    \sum_{i \ne j}^{N_\text{R, T}} \left< {\cal F}_{ij}^{(m)} (0) {\cal F}_{ij}^{(m)} (t)  \right>,

where :math:`N_\text{R}` and :math:`N_\text{T}` are partial ensembles,
where R denotes the rotational and T the translational relaxation modes,
see Ref :cite:`singerMolecularDynamicsSimulations2017` for more details.

.. bibliography::
   :style: unsrt