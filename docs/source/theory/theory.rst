
Theory
======

Dipolar relaxation
------------------

.. container:: justify

    The system of interest here is an ensemble of identical spins characterized by a gyromagnetic ratio  :math:`\gamma_I`
    and spin quantum number :math:`I`. For ¹H, the most abundant isotope of hydrogen, :math:`I = 1/2` and :math:`\gamma_I = 26.752` rad/T/s.
    For ¹³C, a natural and stable isotope of carbon, :math:`I = 1/2` and :math:`\gamma_I = 6.728` rad/T/s :cite:`kowalewskiNuclearSpinRelaxation2006`.

    One assumption is that the cross-correlation terms can be neglected, see Ref :cite:`lippensRelaxationTimeWater1993`.

    When the spin-lattice relaxation is dominated by fluctuations of the magnetic dipole-dipole interaction,
    as is the case for protons in molecular systems, the rates :math:`R_1 (\omega)` and :math:`R_2 (\omega)` are
    related to the spectral densities :math:`J(m)(\omega)` of these fluctuations via the
    Bloembergen-Purcell-Pound (BPP) equations :cite:`bloembergenRelaxationEffectsNuclear1948`:

.. math::
    :label: eq_BPP

    R_1 (\omega) & = & K \left[ J^{(1)} (\omega) + J^{(2)} (2 \omega) \right],

    R_2 (\omega) & = & K \left[ J^{(0)} (0) + 10 J^{(1)} (\omega) + J^{(2)} (2 \omega) \right] / 4,

where

.. math::

    K = \dfrac{3}{2}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

.. container:: justify

    where :math:`\mu_0` is the vacuum permeability, and :math:`\hbar` the Planck constant (divided by :math:`2 \pi`).
    The constant :math:`K` has the units of :math:`\text{m}^6/\text{s}^2`.
    The spectral densities :math:`J^{(m)} (\omega)` in Eq. :eq:`eq_BPP` can be obtained as the Fourier transform
    of the autocorrelation functions :math:`G^{(m)}(\tau)`

.. math::

    J^{(m)} (\omega) = \int_0^\infty G^{(m)} (\tau) \cos(\omega \tau) \mathrm d \tau.

.. container:: justify

    The spectral densities are a measure of the distribution of the fluctuations of :math:`G^{(m)}(\tau)`
    among different frequencies, so they provide information on the distribution of the power available
    for causing spin transitions among different frequencies.
    The autocorrelation functions :math:`G^{(m)}(\tau)` read

.. math::

    G^{(m)} (\tau) = \left< F_2^{(m)} [\textbf{r}_{ij} (t)] F_2^{*(m)} [\textbf{r}_{ij} (0)] \right>

.. container:: justify

    where :math:`F_2^{(m)}` are some complex functions of the vector :math:`\textbf{r}_{ij}` between the spin pairs,
    with norm :math:`r_{ij}` and orientation :math:`\Omega_{ij}` with respect to a reference applied magnetic
    field, assumed to be in the :math:`\textbf{e}_z` direction. The functions :math:`F_2^{(m)}` read 

.. math::

    F_2^{(m)} [\textbf{r}_{ij} (t)] = \alpha_m \dfrac{Y_2^{(m)} [\Omega_{ij} (t)]}{r_{ij}^3 (t)}

.. container:: justify

    where :math:`Y_2^{(m)}` are normalized spherical harmonics and
    :math:`\alpha_0^2 = 16 \pi /5`, :math:`\alpha_1^2 = 8 \pi /15`, and :math:`\alpha_2^2 = 32 \pi / 15`.
    Therefore, one can write for the correlation functions:

.. math::

    G^{(m)} (\tau) = \dfrac{\alpha_m^2}{N}
    \sum_i \sum_{j \ne i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)} \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)},

.. container:: justify

    where :math:`N` is the number of spins.

Intra/inter contributions
-------------------------

.. container:: justify

    Intra-molecular and inter-molecular contributions to :math:`R_1` and :math:`R_2`
    can be extracted separately, by splitting the correlation functions as:

.. math::
    :label: G_intra

    G^{(m)}_\text{intra} (t) = \dfrac{\alpha_m^2}{N}
    \sum_i \sum_{j \in M_i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)}
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)},

.. math::
    :label: G_inter

    G^{(m)}_\text{inter} (t) = \dfrac{\alpha_m^2}{N}
    \sum_i \sum_{j \notin M_i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij}^3 (0)}
    \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij}^3 (\tau)},

.. container:: justify

    where :math:`j \in M_i` and  :math:`j \notin M_i` refer to summation on spin from the 
    same molecule as :math:`i`, and from different molecules as :math:`i`, respectively.

    Intra-molecular relaxation is usually attributed to the rotational motion of the 
    molecule, and inter-molecular relaxation to the translational motion. This is only an 
    assumption that simplify the interpretation of results, and can lead to 
    error :cite:`hubbardTheoryNuclearMagnetic1963`.

Isotropic system
----------------

.. container:: justify

    For isotropic system, the correlation functions are proportional to each others, 
    and only :math:`G^{(0)} (t)` needs to be calculated.

    In that case, :math:`G^{(0)} = 6 G^{(1)}`, and :math:`G^{(0)} = 6 / 4 G^{(2)}`
    :cite:`becherMolecularDynamicsSimulations2021`, which can easily be checked, for instance
    for a bulk water system:

.. image:: ../figures/best-practices/proportionality-dark.png
    :class: only-dark
    :alt: NMR results obtained from the LAMMPS simulation of water

.. image:: ../figures/best-practices/proportionality-light.png
    :class: only-light
    :alt: NMR results obtained from the LAMMPS simulation of water

.. container:: justify

    Therefore, the spectrums can be calculated as:

.. math::

    R_1 &=&  \frac{K}{6} \left[ J^{(0)} (\omega_0) + 4 J^{(0)} (2 \omega_0) \right],

    R_2 &=& \frac{K}{6} \left[ J^{(0)} (0) + \frac{5}{2} J^{(0)} (\omega_0) + J^{(0)} (2 \omega_0) \right],

.. container:: justify

    which require less computational time and less memory to achieve, as only 

.. math::
    :label: F_2_0

    F_2^{(0)} [\textbf{r}_{ij} (t)] & = & \alpha_m \dfrac{Y_2^{(0)} [\Omega_{ij} (t)]}{r_{ij}^3 (t)}

    & = & \dfrac{3 \cos^2 \theta_\text{ij} (t) - 1}{r_{ij}^3 (t)}

.. container:: justify

    needs to be evaluated.

