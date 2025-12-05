\page Examples_channeling Category "channeling"

 Examples in this directory are dedicated to various coherent effects
 in oriented crystals, in particular, channeling, channeling radiation,
 coherent bremsstrahlung, coherent pair production etc. as well as their
 various applications.

\ref Examplech0

 This example shows how channeling in bent crystal can be simulated
 in Geant4 using G4Channeling process. The example simulates the channeling 
 of 400 GeV/c protons in bent Si crystal. It has been moved into channeling/ch0 
 from the channeling folder.

\ref Examplech1

 This example is an easy demonstration of the minimum requirements necessary
 to integrate the G4ChannelingFastSimModel and the G4BaierKatkov model 
 into a project in order to simulate the physics of channeling and 
 channeling radiation/coherent bremsstrahlung.

\ref Examplech2

 This example is an enhanced version of ch1, providing the user with 
 the full functionality of both the G4ChannelingFastSimModel and G4BaierKatkov, 
 with parameters set up via a macro, in order to simulate the physics of 
 channeling and channeling radiation/coherent bremsstrahlung, and enhanced output.
 The example can be exploited for a wide range of cases to study coherent effects in
 a straight, bent or periodically bent crystal (crystalline undulator).
 
\ref Examplech3

 This example is an easy demonstration of the minimum requirements necessary
 to integrate the G4CoherentPairProduction process along with G4ChannelingFastSimModel
 and G4BaierKatkov into a project in order to simulate the physics of electromagnetic
 shower in an oriented crystal. The simulation includes the physics of channeling,
 channeling radiation/coherent bremsstrahlung and coherent pair production.
 The structure of this example is based on ch1, but with the G4CoherentPairProductionPhysics
 process included, as well as different output, a different crystal material, alignment and
 geometry parameters, and a photon beam as the incoming source instead of charged particles.

 \ref Examplech5

 Example ch5 is an application for simulating a positron source.
 Although the conventional approach based on an amorphous target is possible,
 the application is primarily designed to simulate positrons based on oriented crystals.
 In the latter case, both the single-crystal and the hybrid scheme can be investigated.
