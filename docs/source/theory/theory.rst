
Theory
======

Single spin properties
----------------------

NMR makes use of the behavior of quantum mechanical spins in an external magnetic field :cite:`kowalewskiNuclearSpinRelaxation2006`. 
Every nuclear species is characterized by a nuclear spin quantum number, :math:`I`, which can attain integer
or half-integer values (A nucleus with :math:`I = 0` has no spin angular momentum and is not active in NMR.).
The nuclear spin quantum number is related to the eigenvalue of
the total spin angular momentum operator, :math:`\hat I^2`, through:

.. math::

    \hat I^2 \Psi = I (I + 1 ) \Psi

where :math:`\Psi`` is an eigenfunction of the operator :math:`\hat I^2`. Another important operator is the 
:math:`z` component of the spin angular momentum vector, 

.. math::

    \hat I_z \Psi = m \Psi

where the quantum number :math:`m` can attain values ranging between :math:`-I` and :math:`I`, in steps of one.
For spin with :math:`I = 1/2`, :math:`m` has two possible eigenvalues, `-1/2` and `1/2`. The spin angular momentum
is related to the nuclear magnetic moment :math:`\mu_z`, 

.. math::

    \mu_z = \gamma_I I_z

where :math:`\gamma_I` is gyromagnetic ratio (for ¹H, :math:`I = 1/2` and :math:`\gamma_I = 26.752` rad/T/s,
for ¹³C, :math:`I = 1/2` and :math:`\gamma_I = 6.728` rad/T/s :cite:`kowalewskiNuclearSpinRelaxation2006`).

A nuclear spin and a magnetic field :math:`B` of magnitude :math:`B_0` interact through
the so-called Zeeman interaction, which is described by the Zeeman Hamiltonian (or the total energy operator):

.. math::

    \hat H_z = - \gamma_I B_0 \hat I_z

For :math:`I = 1/2`, the Zeeman Hamiltonian has two eigenvalues, :math:`E_{ \pm 1/2} = \mp \frac{1}{2} \gamma_I B_0`
(in angular frequency units). The difference between the two eigenvalues, :math:`\omega_0 = \gamma_I B_0`,
is called the Larmor frequency. Note that quantum mechanics does not require the system to be in a specific
eigenstate of the Hamiltonian, and the system can also exist in a superposition state.

Crucially for NMR relaxation, the time evolution of the wave function for a quantum system is given by
the time-dependent Schrödinger equation:

.. math::

    \partial_t \Psi (t) = - i \hat H \Psi (t)

Spin relaxation
---------------

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

    G^{(m)} (t) = K \dfrac{1}{N}
    \sum_{i \ne j}^{N} \left< {\cal F}_{ij}^{(m)} (0) {\cal F}_{ij}^{(m)} (t)  \right>,

where :math:`N` is the number of spin pairs, and

.. math::

    K = \dfrac{3}{2}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

where :math:`\mu_0` is the vacuum permeability and :math:`I = 1/2` the
spin quantum number. The constant :math:`K` has the units of :math:`\text{m}^6/\text{s}^2`, and therefore 
the functions :math:`G^{(m)}` the units of :math:`\text{s}^{-2}`. Spectral densities are obtained from the
Fourier transforms of the correlation functions, 

.. math::

    J^{(m)} (\omega) = \int_\infty^\infty G^{(m)} (t) \mathrm e^{- i \omega t} \mathrm dt 

from which the relaxation rates can be calculated as

.. math::

    R_1 &=&  J^{(1)} (\omega_0) + J^{(2)} (2 \omega_0),

    R_2 &=& \dfrac{1}{4} \left[ J^{(0)} (0) + 10 J^{(1)} (\omega_0) + J^{(2)} (2 \omega_0) \right],

were :math:`\omega_0 = \gamma B_0` is the Larmor frequency.

Intra/inter contributions
-------------------------

Intra-molecular and inter-molecular contributions to :math:`R_1`
can be extracted separately, by splitting the correlation functions as:

.. math::

    G^{(m)}_\text{R, T} (t) = K \dfrac{1}{N_\text{R, T}}
    \sum_{i \ne j}^{N_\text{R, T}} \left< {\cal F}_{ij}^{(m)} (0) {\cal F}_{ij}^{(m)} (t)  \right>,

where :math:`N_\text{R}` and :math:`N_\text{T}` are partial ensembles,
where R denotes the rotational and T the translational relaxation modes,
see Ref :cite:`singerMolecularDynamicsSimulations2017` for more details.

Isotropic system
----------------

For isotropic system, the correlation functions are proportional to each others, 
and only :math:`G^{(0)} (t)` needs to be calculated.

In that case, :math:`G^{(0)} = 6 G^{(1)}`, and :math:`G^{(0)} = \frac{3}{2} G^{(2)}`.

For isotropic system, spectrums can be calculated as:

.. math::

    R_1 &=&  \frac{1}{6} \left[ J^{(0)} (\omega_0) + 4 J^{(0)} (2 \omega_0) \right],

    R_2 &=& \frac{1}{6} \left[ J^{(0)} (0) + \frac{5}{2} J^{(0)} (\omega_0) + J^{(0)} (2 \omega_0) \right].

.. bibliography::
   :style: unsrt