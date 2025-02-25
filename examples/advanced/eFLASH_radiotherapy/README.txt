
 	     =========================================================
                     Text version of the eFLASH_radiotherapy README file
             =========================================================

Authors: 
Jake Harold Pensavalle (1,2), Said Ahmad (3),  Giuliana Milluzzo (3) and Francesco Romano (3)
 

 (1) Fisica Sanitaria, Azienda Ospedaliero Universitaria Pisa AOUP, ed.18 via Roma 67 Pisa, Italy

 (2) Istituto Nazionale di Fisica Nucleare, Sezione di Pisa, Largo B. Pontecorvo 3 I Pisa, Italy

 (3) Istituto Nazionale di Fisica Nucleare, Sezione di Catania, Via Santa Sofia 64 Catania, Italy


---Contacts:

Jake Pensavalle: jake.pensavalle@pi.infn.it

Giuliana Milluzzo: giuliana.milluzzo@ct.infn.it

Francesco Romano: francesco.romano@ct.infn.it


WHAT IT IS, WHAT IT DOES AND WHAT IT WILL PROVIDE

This Example "FLASH" simulates the beamline and energy spectra based on the Triode Electron Gun Equipped ElectronFlash Manufactured by Sordina Iort Technologies S.p.A. available at the Centro Pisano Flash Radiotherapy (CPFR) at Pisa.  The linac is able to provide low-energy electron flash beams and is currently the first machine installed in Italy capable of reaching  FLASH dose-rate (MGy/s) and extremely high dose per pulse DPP). 


=====================BEAMLINE===============================
The beamline (set along the x-axis) is modelled in the Applicator class, starting from the titanium windows (which is after the RF cavities not modelled in this simulation) ending at the PMMA applicator which acts as a passive collimation of the electron beam.
The applicator has 5mm thick walls and the user can set the inner and outer radius and the total length of the applicatior in the Applicator.cc class by modifying the appropriate parameters (OuterRadiusFirstApplicatorFlash and hightFinalApplicatorFlash)  by means of the SetApplicatorLength and SetOuterRadius functions. The default is the 100mm diameter applicator (inner diameter), but one can 
simulate the others as follows:

1) Applicator (internal diameter) 30mm : OuterRadiusFirstApplicatorFlash = 17.5 mm,  hightFinalApplicatorFlash = 175 mm

2) Applicator (internal diameter) 40mm : OuterRadiusFirstApplicatorFlash = 22.5 mm,  hightFinalApplicatorFlash = 200 mm

2) Applicator (internal diameter) 50mm : OuterRadiusFirstApplicatorFlash = 27.5 mm,  hightFinalApplicatorFlash = 225 mm

Optionally one can set the length of the applicator to 0 and simulate the beam at  the exit window, an additional 3cm of air are between the phantom to account for the applicator inserter hooks on the machine. 

=====================DETECTOR CONSTRUCTOR===============================

The FlashDetectorConstruction Class creates a cubic water phantom (default dimension 30x30x30 cm3) and a detector volume. 
The phantom is placed in contact with the applicator as a typical experimental setup as default. The Phantom sizes can be easily modified by using the SetPhantomSize function at the beginning of the class.
The user can also add an air gap between the end of the final applicator and the entrance of the water phantom by modifying the parameter (AirGap) and the function SetAirGap (default AirGap=0 cm). By default the entrance surface is set to the coordinates (0,0,0) as respect to the center of the world. 
The defaul material of the phantom is water (G4_WATER) but the material can be easily changed by using the following macro commands:
Es:
/changePhantom/material G4_AIR
/chagePhantom/update

The user may also select the maximum step length in this class by changing the "maxStep" parameter.

An array of Silicon Carbide (SiC) detectors is also implemented as geometry (default number of detectors=40, default center to center distance=3 mm, default single SiC detector size: 2x2 mm2 area, 10 um thickness for the active volume, 370 um thickness for the substrate layer) inside the water phantom and it is placed at 13 mm from the surface of the water phantom (depth of the maximum dose for 9 MeV electron beams).
The width and the thickness of the active layer and of the substrate can be easily modified by means of the SetDetectorThickness, SetDetector_subThickness, SetDetectorWidth functions at the beginning of the class.
The user can change the number of detectors and center to center distance by modifying the parameters (nDet and fDet_ctc respectively) at the beginning of the class. nDet can be set to 1 to have a single SiC detector.
The SiC array is not activated by default and can be activated through the following macro command:

/changeDetector/activateArray true 

This macro command must be applied before initialization (/run/initialize).

The user may also change the material of the detector trough macro commands: 

/changeDetector/material G4_C
/chageDetector/update


=====================PRIMARY GENERATOR=====================

The source is modelled in the FlashPrimaryGeneratorAction class and fully uses the General Particle Source class (both in .mac commands and directly in c++ code). 
The beam spot is gaussian with a slight angular divergence (1.5 deg) and the electron energy spectrum can be loaded from macro commands  chosing the 9MeV mode (using "/control/execute 9MeVEF.mac") 
or the 7MeV mode (using "/control/execute 7MeVEF.mac). The geometry of the beam can be set by the user by changing the position, sigma and theta values in 
FlashPrimaryGeneratorAction.cc. 


=====================SCORING MESH=====================

Two scoring meshes are implemented as default in .mac file commands.
The first one (boxMesh_1) is a 120x4x4 mm3 scorer with 120x1x1 number nbin while second scorer (boxMesh_1) is a 4x4x150 mm3 scorer with 1x1x150 number of beam, simulating a similar volume as the Advanced Markus Chamber created in the position of the water phantom to register the dose in Gy delivered in each voxel. Results are printed in the file dose_longitudinal.out and dose_transversal.out at the end of the simulation and represent the depth dose and transversal dose distributions respectively of the electron beams (9 MeV as default) as reported in Di Martino et al., Front. Phys. 11:1268310.
User may change the dimensions and bin number (voxel size) to appropiately score the volume of interest and optimize the simulation (a very finely voxellized volume may increase the simulation time significantly).


=====================SENSITIVE DETECTOR=====================
The FlashSensitiveDetector class is also imlemented for the Silicon carbide detector array. It is only activated if the array is activated through the following command: /changeDetector/activateArray true. 
The output file is .csv file by default. The user may change the output file format (e.g. .root) by modifying the argument of fAnalysisManager->OpenFile inside FlashRunAction::BeginOfRunAction.

An output csv file (output_0_nt_fSensitiveDetector_t*.csv)  will be produced for each thread (in case of multhidreading=ON) and contains the following information for each Silicon carbide detector composing the array and for each step within the sensitive volume:

detector position x, detector position y, detector position z, dose deposited, energy deposited, eventID, parentID, particle name 



=====================ENERGY SPECTRUM AT THE WATER PHANTOM ENTRANCE=====================
The primary particle energy spectrum at the entrance of the phantom (downstream the beamline) is also registered using the FlashSteppingAction class. An ASCII file is automatically written with the following information:

parentid eventid kineticEnergy pos_x pos_y pos_z momentum  cos_x cos_y cos_z

If user is not interested on the retrieval of the energy spectrum can easily comment the line in the SteppingAction class. 


=====================PHYSICS PROCESSES AND PHYSICS LIST =====================

EM Standard option 4 is activated. The user can change the physics list in the physics list class. A production cut for gammas e+ and e- can be applied provided the specification of the Region name defined 
in the FlashDetectorConstruction Class (Default region is the water phantom). 



