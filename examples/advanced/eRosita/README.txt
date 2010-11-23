$Id: README.txt,v 1.4 2010-11-23 20:09:32 pia Exp $
-------------------------------------------------------------------

     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                       eROSITA fluorescence
                       --------------------

Authors:

Dieter Schlosser (pnSensor, Munich), Georg Weidenspointner (MPE Garching
and MPI Halbleiterlabor, Munich), and Maria Grazia Pia (INFN Genova)

References:

M.G. Pia et al., 2009, "PIXE Simulation With Geant4", 
IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649

N. Meidinger et al., 2010, "Development of the focal plane PNCCD
camera system for the X-ray space telescope eROSITA", Nuclear
Instruments and Methods in Physics Research A 624, 321-329

Web sites:

????

Contact persons:

Dieter Schlosser     : dis@hll.mpg.de
Georg Weidenspointner: ggw@hll.mpg.de
Maria Grazia Pia     : Maria.Grazia.Pia@cern.ch


Example description:

This is example is based on simulations of the instrumental background
of the eROSITA X-ray telescope, in particular the strength of
fluorescence lines inside the passive graded Z shield. The set-up
considered in this example consists of a Cu block that is irradiated
with protons. Impact ionization of Cu atoms generates vacancies in
atomic shells. These are then filled by atomic de-excitation,
resulting in the emission of fluorescence photons and Auger electrons
- the PIXE (particle induced X-ray emission) process. In particular
fluorescence photons are then detected with a Si CCD positioned next to
the Cu block. 

The simulated set-up is defined in eRositaDetectorConstruction.
The Cu block is located at position x = y = z = 0 cm. Its dimensions
in x, y, and z are 0.5 cm * 0.5 cm * 3 cm. The CCD, represented by a
slab of Si, is positioned at x = z = 0 cm and y = 2 cm. The CCD dimensions
in x, y, z and are 4 cm * 450 mu_m * 4 cm.

The vertex and initial momenta of the protons are defined in
eRositaPrimaryGeneratorAction. The vertex is at x = 0 cm, y = 2.25
cm, and z = 4 cm. The initial direction of the protons is given by the
vector (0.0, -0.5, -1.0). The initial kinetic energy of the protons is
150 MeV.

The physics processes relevant to this simulation are defined in 
eRositaPhysicsList. The key process for generation of PIXE is 
G4hImpactIonisation.

Instructions on how to build and run the example:

- To compile the example:
  % cd eRosita
  % gmake
  If the environment variable G4WORKDIR and has been defined, an executable
  named eRosita will be generated in $G4WORKDIR/bin/$G4SYSTEM

- To run the example interactively:
  Do not forget to define the G4PIIDATA environment variable as appropriate
  to get access to the PIXE data library 
  
  To run without visualisation:
  $G4WORKDIR/bin/$G4SYSTEM/eRosita eRosita.in
  
  To run with visualisation:
  Copy eRosita/vis.mac to $G4WORKDIR/bin/$G4SYSTEM
  Then go to $G4WORKDIR/bin/$G4SYSTEM and run eRosita
  The main is written such that the example will automatically read 
  the macro vis.mac when running in interactive mode.


Description of output:


