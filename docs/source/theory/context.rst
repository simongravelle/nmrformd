
Context
=======

The measurement of NMR relaxation quantities allow for detailed studies of molecular motions
on time scales ranging from microseconds to minutes in systems as diverse as gases,
liquids, gels, polymers, adsorbed liquids, or solids
:cite:`goreNMRRelaxationWater1989, greiner-schmidSelfDiffusionCompressed1991`,
as well as proteins and other biological systems
:cite:`jacobsonProtonMagneticResonance1954a, rorschachProteinDynamicsNMR1986`.

The key ingredient needed for an accurate description of the nuclear spin relaxation
of 1H in soft matter systems is the accurate description of the random rotational and
translational motions that molecules undergo, which makes classical molecular dynamics (MD)
simulations a natural choice.  For instance NMR relaxation have been used to characterize Lennard-Jones fluid 
MD simulations :cite:`odeliusIntermolecularDipoleDipoleRelaxation1993, grivetNMRRelaxationParameters2005`,
water and other small molecules :cite:`lippensRelaxationTimeWater1993, calero1HNuclearSpin2015, singerMolecularDynamicsSimulations2017, philipsProtonNMRRelaxation2019, singerElucidatingNMRRelaxation2020`.
MD simulations also allow the study of molecules confined within
nanoporous materials :cite:`khudozhitkovDynamicsPropenePropane2020, gravelleNMRInvestigationWater2023a`


In addition to classical MD, Ab initio MD has also been used to extract NMR relaxation time
from water  :cite:`calero1HNuclearSpin2015`. Ab initio and its variants are 
particularly useful in the case where the relaxation mechanisms involve quadrupolar interactions
:cite:`philipsQuadrupolarNMRRelaxation2020, chubakNMRRelaxationRates2021`,
which is beyond the scope of the present contribution.
Monte carlo simulations have also been used :cite:`friesMonteCarloCalculation1983`,
although in that case some care must be taken to extract time dependant quantities
such as autocorrelation functions :cite:`huitemaCanMonteCarlo1999`.  Recently,
the possibility to calculate NMR relaxation rate from coarse grained models was also
discussed :cite:`gravelleAssessingValidityNMR2023`. 

Authors have also directly combined simulations and experiments, for instance to study
glass transition phenomenon of glycerol :cite:`becherMolecularDynamicsSimulations2021`,
water confined within salt crusts :cite:`gravelleNMRInvestigationWater2023a`
