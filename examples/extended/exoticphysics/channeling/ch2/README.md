\page Examplech2 Example ch2

\author Alexei Sytov, Gianfranco Paternò - INFN Ferrara Division (Italy) \n
 sytov@fe.infn.it, paterno@fe.infn.it

## INTRODUCTION
Example ch2 is an enhanced version of ch1, providing the user with the full functionality of
both the G4ChannelingFastSimModel and G4BaierKatkov, with parameters set up via a macro,
in order to simulate the physics of channeling and channeling radiation/coherent bremsstrahlung.

The example can be exploited for a wide range of cases to study coherent effects in
a straight, bent or periodically bent crystal (crystalline undulator). Channeling
physics in ch2 is active for protons, ions, muons, pions, electrons and their antiparticles.
Any other charged particle can also be activated.

The example contains also other setups for specific applications.

## DESCRIPTION
The setup of the example ch2 in run.mac is identical to ch1. As ch1, this example includes a bent crystal
and a detector positioned behind it. Like ch1, it is based on the experiments on
channeling [1] and channeling radiation [2] in a bent crystal, carried out at
Mainz Mikrotron MAMI with 855 MeV electrons. The experimental validation of
G4ChannelingFastSimModel is described in [3].

However, since ch2 parameters are fully set up in the macro run.mac, this example
is quite flexible and can be easily adapted for entirely different cases.

In addition more specific macros were created to supply users with the setups related to the applications.
These macros partially exploit the model defaults to simplify the example. They include:
-# run_Bent_Crystal_Deflection_Radiation.mac - reduced version (some commands setting defaults deleted) of run.mac but with an identical setup.
-# run_Bent_Crystal_HE_Deflection.mac - an example of particle deflection in a bent crystal at high energies.
-# run_Positron_Source.mac - a simplified example of a positron source within a single W target.
-# run_Radiation.mac - an example of a radiation source in a straight crystal.

A description of all the available options is provided in run.mac and partially in other macros.
It includes crystal and detector geometry, activation flags for
G4ChannelingFastSimModel and G4BaierKatkov and various options. 

The example also provides detailed descriptions of various options for 
G4ChannelingFastSimModel and G4BaierKatkov, which can adjust model parameters 
depending on the specific case (see DetectorConstruction::ConstructSDandField()).

The front surface of the crystal is placed at z=0 (with z as the beam direction), 
while the front position of the detector can be set up via run.mac.

The output is recorded into the file results.root as a set of root ntuples.
These ntuples include:
-# crystal: particles recorded at the crystal entrance,
-# detector_primaries: primaries recorded at the detector entrance AND passed through the crystal.
-# detector_photons: photons recorded at the detector entrance produced by primaries passed through the crystal.
-# detector_sedondaries: secondaries recorded at the detector entrance produced by primaries passed through the crystal.
-# missed_crystal: all the particles missed the crystal, however, entering the detector, if any.

The format of every ntuple includes the following 10 variables (columns):

- "eventID", "volume", "x", "y", "angle_x", "angle_y", "Ekin", "particle", "particleID", "parentID"

The variables represent:
-# the event number within the run (column 0),
-# the volume, either the crystal or the detector (column 1),
-# the coordinate (x,y) and the angles (x'=dx/dz, y'=dy/dz) of the impinging particles (columns 2-6),
-# the kinetic energy of the particle (column 7),
-# the particle name (column 8),
-# the particle ID (column 9),
-# the parent ID of the particle (column 10).

For convenience for detector_primaries were added four more variables:
-# the incoming angle x at the crystal entrance,
-# the deflection angle x (the difference between the angle at the detector and the incoming angle),
-# the incoming angle y at the crystal entrance,
-# the deflection angle y (the difference between the angle at the detector and the incoming angle).

These four variables are especially useful for the studies of deflection of primary particles.

To visualize these data, one should either open results.root using root TBrowser or use the python script analysis_ch2.py or its identical version in the jupyter notebook format analysis_ch2.ipynb.

The output data also includes the spectrum of photons using the data produced inside the Baier-Katkov method.
This spectrum requires nearly 2 order on magnitude less data, then the collection of gamma produced in Geant4 as secondaries.
It is very useful especially if the goal is to produce only the spectrum of radiation. This spectrum is normalized on the
total radiation emission probability, which is an equivalent to 1/Nprimaries dN_photon/dE_photon normalization.

Moreover, it is possible to set up a round virtual collimator - an angular selection of photons in the Baier-Katkov method.
This is extremely useful for coherent bremsstrahlung simulation.

CAUTION: though the Baier-Katkov spectrum should identically coincide with the spectrum by secondary photons, sometimes
it may be less accurate, since it is updated only after every setNSmallTrajectorySteps (see run.mac). Moreover,
the virtual collimator does not take into account the transverse positions of particles. Therefore, it is recommended
to use the Baier-Katkov spectrum at low statistics for preliminary researches and optimization while the secondaries produced at
high statistics as a final result.

CAUTION: the angular center of virtual collimator coincides with the global z direction.

The spectrum is produced as a text file, containing the photon energies in the first column and the corresponding spectrum value in the second one.
Note, the bins are not equidistant, they are sampled according to the bremsstrahlung spectrum, with the bin size proportional to 1/E_photon.

## REFERENCES

-# A. Mazzolari et al. <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.112.135503">Phys. Rev. Lett. 112, 135503 (2014).</a>
-# L. Bandiera et al. <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.025504">Phys. Rev. Lett. 115, 025504 (2015).</a>
-# A. Sytov et al. <a href="https://link.springer.com/article/10.1007/s40042-023-00834-6"> Journal of the Korean Physical Society 83, 132–139 (2023).</a>
