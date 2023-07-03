
Theory (UNDER CONSTRUCTION)
===========================

When the SLR is dominated by fluctuations of the magnetic dipole-dipole interaction,
as is the case for protons in molecular systems, the rate :math:`R_1 (\omega)` is
related to the spectral densities :math:`J(m)(\omega)` of these fluctuations via the
Bloembergen-Purcell-Pound (BPP) equation :cite:`bloembergenRelaxationEffectsNuclear1948`:

.. math::
    :label: eq_BPP

    R_1 (\omega) = K \left[ J^{(1)} (\omega) + J^{(2)} (2 \omega) \right],

where

.. math::

    K = \dfrac{3}{2}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

where :math:`\mu_0` is the vacuum permeability, :math:`\hbar` the Planck constant (divided by :math:`2 \pi`),
:math:`\gamma_I` is the gyromagnetic ratio (for ¹H, :math:`I = 1/2` and :math:`\gamma_I = 26.752` rad/T/s,
for ¹³C, :math:`I = 1/2` and :math:`\gamma_I = 6.728` rad/T/s :cite:`kowalewskiNuclearSpinRelaxation2006`), and
:math:`I = 1/2` the spin quantum number. The constant :math:`K` has the units of :math:`\text{m}^6/\text{s}^2`.

The spectral densities :math:`J^{(m)} (\omega)` in Eq. :eq:`eq_BPP` can be obtained as the Fourier transform
of the autocorrelation functions :math:`G^{(m)}(\tau)`

.. math::

    J^{(m)} (\omega) = \int_0^\infty G^{(m)} (\tau) \cos(\omega \tau) \mathrm d \tau.

The spectral densities are a measure of the distribution of the fluctuations of :math:`G^{(m)}(\tau)`
among different frequencies, so they provide information on the distribution of the power available
for causing spin transitions among different frequencies.
The autocorrelation functions :math:`G^{(m)}(\tau)` read

.. math::

    G^{(m)} (\tau) = \dfrac{\alpha_m^2}{N}
    \sum_i \sum_{j \ne i} \dfrac{Y_2^{(m)} [\Omega_{ij} (0)]}{r_{ij} (0)} \dfrac{Y_2^{*(m)} [\Omega_{ij} (\tau)]}{r_{ij} (\tau)},

where :math:`\alpha_0^2 = 16 \pi /5`, :math:`\alpha_1^2 = 8 \pi /15`, and :math:`\alpha_2^2 = 32\pi /5`,
:math:`N` is the number of spins, and :math:`Y_2^{(m)}` are normalized spherical harmonics. :math:`G^{(m)}(\tau)`
are functions of the vector :math:`\textbf{r}_{ij}` between the spin pairs, which is 
characterised by its norm :math:`r_{ij}` and orientation :math:`\Omega_{ij}`
with respect to a reference applied magnetic field, assumed to be in the :math:`\textbf{e}_z` direction.

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

Magnetization
-------------

The magnitude of the magnetization of an ensemble of
N non-interacting spins at thermal equilibrium reads 

.. math::

    M_0 = \dfrac{N \gamma_I^2 \hbar^2 I (I + 1) B_0}{3 k_B T}.

The evolution of the magnetization vector with time is given by the phenomenological 
Bloch equations:

.. math::

    \dfrac{d M_z}{d t} = \dfrac{M_0 - M_z}{T_1},

    \dfrac{d M_y}{d t} = - M_x \omega_\text{off} \dfrac{M_y}{T_2}

    \dfrac{d M_x}{d t} = M_y \omega_\text{off} \dfrac{M_x}{T_2},

where :math:`\omega_\text{off}` is the frequency offset between the applied radiofrequency
and the Larmor frequency.

Alternative presentation
------------------------

One consider the autocorrelation function of a complex function Y:

.. math::

    G (\tau) = < Y (t) Y^* (t + \tau)>

with :math:`G(0) = \sigma` and 

.. math::

    \lim_{t \to \infty} G (\tau) = 0

Thus, we expect a general time-correlation function to be a decaying function of time,
with an initial value given by the variance of Y. A reasonable choice is:

.. math::

    G(\tau) = G(0) \exp(- | \tau | /\tau_c) 

where :math:`\tau_c` is the correlation time, which is a measure of the time scale of oscillations of the random process
or a measure of the persistence of the correlation between values of :math:`Y(t)` at different points in time.

Spectral density can be obtained as the Fourier transform of :math:`G(t)`:

.. math::

    J (\omega) = 2 \int_0^\infty G(\tau) \exp(- i \omega \tau) \mathrm d \tau.

The spectral density is a measure of the distribution of the fluctuations of :math:`Y(t)` among different frequencies,
so they provide information on the distribution of the power available for causing spin transitions among different frequencies.
The spectral density of an exponentially decaying correlation function is a Lorentzian:

.. math::

    J (\omega) = G(0) \dfrac{2 \tau_c}{1 + \omega^2 \tau_c^2}

In most cases, the spectral densities are linear combinations of Lorentzian functions.

The fundamental molecular dynamic quantities of primary interest for NMR are time-correlation
functions for rank-2 spherical harmonics of the pair of angles specifying the direction of a
given molecule-fixed axis with respect to the laboratory frame.

The functions :math:`Y` reads :cite:`bloembergenRelaxationEffectsNuclear1948`

.. math::

    Y_{0j} & = & \dfrac{1 - 3 \cos^2 \theta_\text{ij} }{r_{ij}^3}

    Y_{1j} & = & \dfrac{ \sin \theta_\text{ij} \cos \theta_\text{ij} \exp{i \phi_{ij}} }{r_{ij}^3}

    Y_{2j} & = &  \dfrac{ \sin^2 \theta_\text{ij} \exp{2 i \phi_{ij}} }{r_{ij}^3}

Dipolar relaxation
------------------

Assuming that two nuclear magnetic moments or magnetic dipoles, :math:`\mu_1` and :math:`\mu_2` are close in space.
The field created by the dipole :math:`\mu_2` reads

.. math::

    \textbf{B}_\text{loc} (\mu_2) = - \dfrac{\mu_0}{4 \pi r^3} \left( \mu_2 - 3 \dfrac{\textbf{rr}}{r^2} \cdot \mu_2 \right)

where :math:`\mu_0` is the permeability of vacuum, :math:`r` is the distance from the origin and :math:`\textbf{rr}` a tensor. 
The classical dipole-dipole energy is 

.. math::

    \textbf{E}_\text{DD} = \dfrac{\mu_0}{4 \pi r^3} \left( \mu_1 \cdot \mu_2 - 3 \mu_1 \cdot \dfrac{\textbf{rr}}{r^2} \cdot \mu_2 \right)

Here :math:`\textbf{r}` is the vector connecting the two dipoles. The quantum mechanical counterpart is 
obtained by replacing the magnetic dipoles by :math:`\gamma_I \hbar \hat{\textbf{I}}` and :math:`\gamma_S \hbar \hat{\textbf{S}}`,

.. math::

    \hat{\textbf{H}}_\text{DD} = - \dfrac{\mu_0 \gamma_I \gamma_S \hbar}{4 \pi r^3} \left( 3 \hat{\textbf{I}} \cdot \dfrac{\textbf{rr}}{r^2} \cdot \hat{\textbf{S}}
    - \hat{\textbf{I}} \cdot \hat{\textbf{S}} \right) = b_\text{IS} \hat{\textbf{I}} \cdot \textbf{D} \cdot \hat{\textbf{S}},

where :math:`b_\text{IS}` is the dipole-dipole coupling constant and :math:`\textbf{D}` is the dipolar tensor, which in spherical polar coordinate reads: 

.. math::

    \textbf{D} = \begin{pmatrix}
            3 \sin^2 \theta \cos^2 \phi - 1 & 3 \sin^2 \theta \cos \phi \sin \phi &  3 \sin \theta \cos \theta \cos \phi  \\
            3 \sin^2 \theta \cos \phi \sin \phi & 3 \sin^2 \theta \sin^2 \phi - 1 & 3 \sin \theta \cos \theta \sin \phi \\
            3 \sin \theta \cos \theta \cos \phi & 3 \sin \theta \cos \theta \sin \phi & 3 \cos^2 \theta - 1
        \end{pmatrix}

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

In that case, :math:`G^{(0)} = 6 G^{(1)}`, and :math:`G^{(0)} = 6 / 4 G^{(2)}` :cite:`becherMolecularDynamicsSimulations2021`.

For isotropic system, spectrums can be calculated as:

.. math::

    R_1 &=&  \frac{1}{6} \left[ J^{(0)} (\omega_0) + 4 J^{(0)} (2 \omega_0) \right],

    R_2 &=& \frac{1}{6} \left[ J^{(0)} (0) + \frac{5}{2} J^{(0)} (\omega_0) + J^{(0)} (2 \omega_0) \right].

The case of small molecules
---------------------------

Small molecules in low-viscosity solutions typically have rotational correlation times of a few tens of
picoseconds or less. In that case the extreme narrowing conditions usually prevail, therefore :math:`J_2(\omega) = J_2(0)`.

.. bibliography::
   :style: unsrt

