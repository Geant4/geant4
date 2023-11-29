--------------------------------------------------------------------------------

     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                      Xray_SiliconPoreOptics
                        ------------------
        	        P.Dondero (1), R.Stanzani (1)
                              Apr 2023

 1. Swhard S.r.l, Genoa (GE), Italy.

--------------------------------------------------------------------------------

 Contacts: paolo.dondero@cern.ch, ronny.stanzani@cern.ch

--------------------------------------------------------------------------------
 Acknowledgements: example developed within the ESA AREMBES Project, Contract n.
 4000116655/16/NL/BW. Valentina Fioretti provided the simplified mass model, as
 described in [1].
--------------------------------------------------------------------------------

 Xray_SiliconPoreOptics is an example of the application of Geant4 in a space
 environment.
 The geometry used in this example represents a single reflective pore used to
 simulate on a smaller scale the effect of the millions of pores forming the
 mirror of the ATHENA Silicon Pore Optics (SPO), as described in [1].
 The main purpose of the simulation is the estimation of the induced residual
 background at the pore exit caused by proton scattering at grazing angles
 (<1deg).
 Reflection steps inside the pore and relevant information are saved on a .root
 file for subsequent analysis [2]. For execution time optimization purposes, only
 particle steps respecting specific conditions (e.g. reflection length and volume
 name) are stored.
 An example of ROOT-based analysis of the output file is included
 ("./analysis/analysis.C") and can be used to obtain basic data representations.
 Xray_SiliconPoreOptics implements a physics list dedicated to space radiation
 interactions, developed within the ESA AREMBES Project for the ATHENA mission,
 called Space Physics List (SPL).
 The example shows a way to optimize the simulation's execution time
 and output size by selectively saving data based on specific combined conditions
 (e.g. position, eventID and process name).
 NOTE: in a multiple-run session, the last run always overrides the root file.

1 - GEOMETRY

 The geometry is given in the GDML format, and consists of a single Silicon pore
 aligned to the ideal optics symmetry axis of the SPO [1], i.e., the Z-axis of
 the Geant4 reference system. The pore has the following parameters:
  - length: ~203.0 mm
  - pore entrance size: ~0.83x0.61 mm
  - pore thickness: 0.17 mm
 Three volumes (DummyEntrance, DummyExit and DummySphere) are used to save the
 state of the particles as they pass.

2 - INPUT FLUX

 100keV protons are emitted with a Cosine-law distribution from a planar surface
 (same dimensions of the pore) at 1mm above the entrance, within a cone of +-1 deg
 aperture, as described in [1].

3 - PHYSICS LIST

 This example implements a dedicated physics list called "Space Physics List",
 developed within the ESA AREMBES Project. This physics list has been designed
 focusing on the ATHENA physics processes, but contains high precision
 models that can be used in a more general space application.
 In details, this physics list provides a custom electromagnetic part combined
 with the QBBC hadronic physics list.
 In addition, the G4EmStandardSS Physics List is used to simulate the single
 scattering inside the pore, as it is associated to a specific region
 from the macro file.
 In general, the use of SS only in selected regions allows the simulation to
 reduce CPU consumption in the majority of the volumes and be very accurate in
 the desired ones.
 The default production cuts are selected for all volumes, i.e. 1mm.

4 - HOW TO RUN THE EXAMPLE

 Compile code and execute Xray_SiliconPoreOptics in 'batch' mode from the macro
 file:
 	./XraySiliconPoreOptics run01.mac
 For this example, the multi-thread (MT) capability of Geant4 is enabled by
 default.
 To specify the desired number of threads, the user can use the command
 "/run/numberOfThreads" in "run01.mac". To show the output from a single thread
 in the terminal, the user can use the "/control/cout/ignoreThreadsExcept
 {THREADNUM}" command.

5 - STEPPING

 Within the "SteppingAction" class relevant information about the particle's
 state are stored in Tuples [2], defined in the "HistoManager" class.
 The tuples contain the following information:
  1. event ID
  2. volume name
  3. track ID
  4. coordinates (x,y,z)
  5. angles (theta, phi)
  6. process name
  7. parent ID
  8. the number of internal reflections whenever the particle reaches one of the
     dummy volumes defined above.

6 - ANALYSIS

 Xray_SiliconPoreOptics provides an analysis macro example (analysis.C) to
 visualize data in the following representations:
  - a histogram for the normalized efficiency for Theta and Phi;
  - a pie chart for the number of reflections inside the pore.
 The normalized efficiency serves to observe the angular distribution of the
 exiting protons, normalized over the total entering particles. A proton is
 selected if it enters the first volume (pore entrance), exits from the second
 empty volume (pore exit) and enters the sphere at the detector side (the
 hemisphere below the pore). No pore interaction is required.
 The pie chart reports the number of reflections with the highest probability.

7 - VISUALISATION

   The visualisation manager is set via the G4VisExecutive class in the main()
   function in xray_SiliconPoreOptics.cc.
   The initialisation of the drawing is done via a set of /vis/ commands in the
   macro vis.mac. This macro is automatically read from the main function when
   the example is used in interactive running mode.

--------------------------------------------------------------------------------

References

 [1] Fioretti V et al. "The Geant4 mass model of the ATHENA Silicon Pore Optics
 and its effect on soft proton scattering", Space Telescopes and Instrumentation
 2018: Ultraviolet to Gamma Ray. Vol. 10699. SPIE, 2018.
 [2] BRUN, René, et al. "The ROOT Users Guide". CERN, http://root.cern.ch, 2003.

