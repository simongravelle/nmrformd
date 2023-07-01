
Isotropic systems
=================

NMR relaxation time :math:`T_1` is calculated from 
the autocorrelation function of fluctuating magnetic 
dipole-dipole interactions (see |mcconnell1987| and |bloembergen1948| for more details),

.. math::

    G_\text{R,T} (t) = \dfrac{1}{N_\text{R,T}} \sum_{i \ne j}^{N_\text{R,T}} 
    \left< {\cal F}_{ij} (t + \tau) {\cal F}_{ij} (\tau)  \right>,

where :math:`\tau` is the lag time. Here :math:`N_\text{R}` and :math:`N_\text{T}`
correspond to partial ensembles for intramolecular and itermolecular interactions,
respectively, where R stands for rotational and T for translational (see |singer17a|). The ensemble average is performed by a 
double summation over spin pair :math:`i` and :math:`j` with :math:`i \ne j`.

The function :math:`{\cal F}_{ij} (t)` reads

.. math::
    
    {\cal F}_{ij} (t) = \sqrt{\dfrac{16}{5 \pi}} \dfrac{Y^0_2 (\theta_{ij} (t))}{r_{ij}^3 (t)},

where :math:`Y^0_2` is the normalised spherical harmonics with :math:`\ell = 2` and :math:`m = 0`,
:math:`r_{ij} (t)` the nuclear spin separation, and :math:`\theta_{ij} (t)` the polar angle
of the direction :math:`\textbf{r}_{ij}` with respect to laboratory axes (assuming that 
the applied static magnetic is parallel to :math:`\textbf{e}_z`).

NMR relaxation times are calculated using

.. math::

    \dfrac{1}{T_{1, \text{R}, \text{T}}} = K \left(J_\text{R,T} (\omega_0) + 4 J_\text{R,T} (2 \omega_0) \right),

were :math:`\omega_0 = \gamma B_0` is the Larmor frequency with :math:`\gamma$` the
gyro-magnetic ratio for :math:`^1`H with spin :math:`I = 1/2`, and 

.. math::

    K = \dfrac{1}{4}\left(\dfrac{\mu_0}{4 \pi}\right)^2 \hbar^2 \gamma^4 I (I+1),

where :math:`\mu_0` is the vacuum permeability.

.. |mcconnell1987| raw:: html

   <a href="https://www.cambridge.org/de/academic/subjects/physics/condensed-matter-physics-nanoscience-and-mesoscopic-physics/theory-nuclear-magnetic-relaxation-liquids?format=PB&isbn=9780521107716" target="_blank">mcconnell1987</a>

.. |bloembergen1948| raw:: html

   <a href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.73.679" target="_blank">bloembergen1948</a>

.. |singer17a| raw:: html

   <a href="https://www.sciencedirect.com/science/article/pii/S1090780717300319" target="_blank">singer17a</a>
