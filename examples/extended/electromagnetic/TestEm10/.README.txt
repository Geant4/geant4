//$Id$

///\file "electromagnetic/TestEm10/.README.txt"
///\brief Example TestEm10 README page

/*! \page ExampleTestEm10 Example TestEm10

 Test for investigation of ionisation in thin absorbers, transition 
 and synchrotron radiations. Default setup for "TestEm10.in" and "TestEm10.large_N.in" is 
 the experiment for XTR with NIM A294 (1990) 465-472 (fig. 11) setup


\section TestEm10_s0 INTRODUCTION
	
  The parameterisations models can be changed simply with:
\verbatim
Idle> /XTRdetector/setModel i (i = 1 to 10)
\endverbatim

 It is NOT needed (and not recommended) to issue the command
 /XTRdetector/update if just the model is changed.

 See macro file "TestEm10.in" for an example.

	
\section TestEm10_s1 GEOMETRY DEFINITION
 
 The "absorber" is a tube made of a given material.                
 
 Three parameters define the absorber :
 - the material of the absorber,
 - the thickness of an absorber,
 - the transverse size of the absorber (the input face is a square). 
    
 The volume "World" contains the "absorber". 
 In this test the parameters of the "World" can be changed , too.

 In addition a transverse uniform magnetic field can be applied.
 
 The default geometry is constructed in Em10DetectorConstruction class,
 but all the parameters can be changed via
 the commands defined in the Em10DetectorMessenger class.
 	
\section TestEm10_s2 AN EVENT : THE PRIMARY GENERATOR
 
 The primary kinematic consists of a single particle which hits the
 absorber perpendicular to the input face. The type of the particle
 and its energy are set in the Em10PrimaryGeneratorAction class, and can
 be changed via the G4 build-in commands of G4ParticleGun class (see
 the macros provided with this example).
 
 A RUN is a set of events.
 	
\section TestEm10_s3 DETECTOR RESPONSE

 Here we test G4PAIionisation , G4IonisationByLogicalVolume and 
 transition radiation processes
 
 A HIT is a record, event per event , of all the 
 informations needed to simulate and analyse the detector response.
 
 In this example a Em10CalorHit is defined as a set of 2 informations:
 - the total energy deposit in the absorber,
 - the total tracklength of all charged particles in the absorber,  
 
 Therefore  the absorber is declared
 'sensitive detector' (SD), see Em10CalorimeterSD, which means they can contribute to the hit.
 
 At the end of a run, from the histogram(s), one can study 
 different physics quantities such as :
                         - angle distribution,
                         - energy deposit,
                         - transmission/backscattering,
                         -  ...
 
 The test contains 10 built-in histograms, which can be activated by
 interactive commands (see the macros runxx.mac for details).

 The histogram files can be viewed using PAW e.g with the commands 
\verbatim
paw> h/file 1 geant4.plot01 or g4.p11
paw> option stat
paw> h/pl 1
\endverbatim
 				
\section TestEm10_s4 PHYSICS DEMO
 
 The particle's type and the physic processes which will be available
 in this example are set in Em10PhysicsList class.
 
 The  messenger classes introduce interactive commands . Using these
 commands the geometry of the detector, the data of the primary
 particle, the limits of the histograms , etc. can be changed. 	

\section TestEm10_s5- HOW TO START ?
 
 - Execute TestEm10 in 'batch' mode from macro files e.g.
\verbatim
% TestEm10   run11.mac
\endverbatim
		
 - Execute TestEm10 in 'interactive' mode with visualization e.g.
\verbatim
% TestEm10
....
Idle> type your commands
....
\endverbatim
		
\section TestEm10_s6 List of the built-in histograms

  - 1.   number of (tracking) steps/event
  - 2.   energy deposit distribution in the absorber (in MeV)
  - 3.   angle distribution of the primary particle at the exit
         of the absorber (deg)
  - 4.   distribution of the lateral displacement at exit(mm)
  - 5.   kinetic energy of the transmitted primaries (MeV)
  - 6.   angle distribution of the backscattered primaries (deg)
  - 7.   kinetic energy of the backscattered primary particles (MeV)
  - 8.   kinetic energy of the charged secondary particles (MeV)
  - 9.   z distribution of the secondary charged vertices (mm)
  - 10.   kinetic energy of the photons escaping the absorber (MeV)


\subsection TestEm10_sub_s61 Using histograms

By default the histograms are not activated. To activate histograms
the environment variable G4ANALYSIS_USE should be defined. For instance
uncomment the flag G4ANALYSIS_USE in GNUmakefile.

To use histograms any of implementations of AIDA interfaces should
be available. See InstallAida.txt
  
*/
