--------------------------------------------------------------------------------

     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                          Xray_TESdetector
                         ------------------
                    P.Dondero (1), R.Stanzani (1)
                              Dec 2022

 1. Swhard S.r.l, Genoa (GE), Italy.

--------------------------------------------------------------------------------

 Contacts: paolo.dondero@cern.ch, ronny.stanzani@cern.ch

--------------------------------------------------------------------------------
 Acknowledgements: example developed within the ESA AREMBES Project,
 Contract n. 4000116655/16/NL/BW. This example is a reduced mass model of the
 Athena X-IFU instrument based on an early configuration which is no more
 applicable for any evaluation. Simone Lotti provided the simplified mass
 model and background derived from those used in [1].
--------------------------------------------------------------------------------

 Xray_TESdetector is an example of the application of Geant4 in a space
 environment. It represents an x-ray detector derived from the X-IFU, the
 X-ray spectrometer designed and developed by the European Space Agency (ESA)
 for use on the ATHENA telescope.
 The detector is a Transition-edge sensor (TES) composed of 317 Bismuth pixels
 arranged in a hexagonal shape and its setup includes different layers of
 shielding, filters and support.
 The main purpose of the simulation is the estimation of the particle radiation
 background impacting on the detector. For execution time optimization purposes,
 only particle steps respecting specific conditions (e.g. hit selected
 volumes close to the detector) are stored on a .root file [2].
 An example of ROOT-based analysis of the output file is included
 ("./analysis/analysis.C") and can be used to obtain basic plots and histograms.
 Xray_TESdetector implements a physics list dedicated to space radiation interactions,
 developed within the ESA AREMBES Project for the ATHENA mission, called Space
 Physics List (SPL).
 Technically, this example shows how to manage a complex geometry obtained
 with advanced detector construction features (e.g., boolean operations,
 parameterisation).
 In addition, the example shows a way to optimize the simulation's execution time
 and output size by selectively saving data based on specific combined conditions
 (e.g. position, eventID and process name).
 NOTE: in a multiple-run session, the last run always overrides the root file.

1 - GEOMETRY DEFINITION

 The geometry consists of a simplified version of the X-IFU detector and is composed of
 the following:
  - the TES array, the backscattering (BSC) and the
 Anti-coincidence detector (ACD);
  - the structural elements supporting the detector (e.g. the cage underneath
 it);
  - the thermal shieldings;
  - the structural elements of the cryostat chamber;
  - a hollow Aluminum sphere schematizing the satellite.
 Detector parameters:
  - Detector thickness: 3 um
  - Number of pixels: 317
  - Detector's shape: regular hexagon
  - Hexagon's apothem: 8.593 mm
 The default geometry is constructed in DetectorConstruction class.
 Alternatively, a GDML file is provided (xray_TESdetector.gdml).
 The position of each pixel is defined by a list of coordinates (x,y)
 contained in "pixelpos.txt".

2 - PARTICLE SOURCE

 The radiation field is composed of galactic cosmic rays (GCR) protons with a
 flux estimated for the L1/L2 Lagrangian points, as described in [1].
 The energies range from 10 MeV to 100 GeV, and the particles are isotropically
 generated on the surface of a sphere surrounding the geometry and randomly
 launched toward its interior. The detector is placed in the center of the
 sphere and the sphere's radius is chosen to avoid intersections with geometry
 elements.

3 - PHYSICS LIST

 This example implements a dedicated physics list called "Space Physics List",
 developed within the ESA AREMBES Project. This physics list has been designed
 focusing on the ATHENA physics processes, but contains high precision
 models that can be used in a more general space application.
 In details, this physics list provides a custom electromagnetic part combined
 with the QBBC hadronic physics list.
 In volumes near the detector, where high precision in the scattering description
 is needed, the use of Single Scattering (SS) model is reccomended, as shown in
 the "run01.mac", through the SetEmRegion command.
 The use of SS only in selected regions allows the simulation to reduce CPU
 consumption in the majority of the volumes and be very accurate near the
 detector.
 The default production cuts are selected for all volumes, i.e. 1mm.

4 - HOW TO RUN THE EXAMPLE

 Compile code and execute the example in 'batch' mode from the macro file:
 	./XrayTESdetector run01.mac
 to launch it with the DetectorConstruction, or:
  	./XrayTESdetector run02.mac
 to launch it by using the provided GDML.
 For this example, the multi-thread (MT) capability of Geant4 is enabled by
 default. To specify the desired number of threads, the user can use the
 command "/run/numberOfThreads" in "run01.mac".

5 - STEPPING

 Within the "SteppingAction" class relevant information about the particle's
 state are stored in Tuples [2], defined in the "HistoManager" class.
 The tuples contain the following information:
  1. event ID
  2. volume name
  3. track ID
  4. coordinates (x,y,z)
  5. angles (theta, phi)
  6. parent ID
  7. pixel number (from the TES array)
  8. step energy deposit
  9. step number
  10. initial kinetic energy
  11. kinetic energy
  12. particle name
  13. pre and post-step names
  14. creator process name

 Tuples are filled with the informations listed above in two cases:
  - when a new particle is generated (both primaries and secondaries);
  - when the particles reach the volumes next to the detector and the
 detector itself.

6 - ANALYSIS

 xray_TESdetector provides an analysis macro example (analysis.C) with several
 predefined histograms:
  - Average energy deposit per pixel (1D);
  - Energy deposit on the detector (2D);
  - Particle count per pixel (1D);
  - Spectra of the primaries on the detector (1D);
  - Total spectra on the detector (1D);
  - distribution of the particles on the detector (1D).

 The first three are used to qualitatively check how the interactions are
 distributed on the detector pixels and what is the average deposit per pixel
 and particle. The 2D histogram for the Energy deposit on the detector shows
 the shape of the detector on the XY plane.
 The spectrum histograms are used to observe the following:
  - the initial energies of the particles (at launch or generation);
  - the energy deposit on the detector;
  - the energy of the step before the impact on the pixel.
 Those information are the starting point to assess the background
 composition and intensity on the detector, and thus optimize the
 detector shielding and background rejection techniques.
 Histograms are managed by the "analysis.C" file.

 7 - VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in xray_TESdetector.cc.
   The initialisation of the drawing is done via a set of /vis/ commands
   in the macro vis.mac. This macro is automatically read from
   the main function when the example is used in interactive running mode.

--------------------------------------------------------------------------------

References

 [1] S. Lotti, S. Molendi, C. Macculi, V. Fioretti, L. Piro et al., "Review of
 the Particle Background of the Athena X-IFU Instrument", The Astrophysical
 Journal, 2021.
 [2] BRUN, René, et al. "The ROOT Users Guide". CERN, http://root.cern.ch, 2003.
