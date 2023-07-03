Best practice (UNDER CONSTRUCTION)
==================================

Box size
--------

Simulation duration
-------------------

Dumping frequency
-----------------

Molecule rigidity
-----------------

Check isotropy
--------------

Although generally true for bulk systems, as showed for glycerol :cite:`becherMolecularDynamicsSimulations2021`,
it can be worth ensuring that the relation

.. math::

    \frac{1}{6} G^{(0)} (t) = G^{(1)} (t) = \frac{1}{4} G^{(2)} (t) 

actually stands. If not, all three correlation functions must be calculated.