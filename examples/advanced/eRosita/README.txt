$Id: README.txt,v 1.6 2010-11-30 13:44:05 pia Exp $
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
G4hImpactIonisation. In the example, proton cross sections based on
the ECPSSR model are used. Cross sections are computed for an energy
range from 1 keV to 200 MeV. The energy thresholds for the production
of fluorescence photons and Auger electrons by PIXE are set to a value
of 250 eV.

The output of the example is an ASCII file named
TrackerPhotonEnergy.out. It contains the energy, in MeV, of every
photon that finds its way to the tracker (or which is created as a
secondary inside the tracker). If the batch mode example is run, the
entries 0.00798467, 0.00800571, and 0.00886534 correspond to Cu
fluorescence lines K_alpha2, K_alpha1, and K_beta1,
respectively. Other photons originate e.g. from bremsstrahlung of
delta rays in Cu. A histogram of the energies in
TrackerPhotonEnergy.out, in particular if generated for more protons,
clearly shows the PIXE photons from Cu on top of a continuous
distribution.


Instructions on how to build and run the example:

- To compile the example:
  % cd eRosita
  % gmake
  If the environment variable G4WORKDIR and has been defined, an executable
  named eRosita will be generated in $G4WORKDIR/bin/$G4SYSTEM

- To run the example:

  Do not forget to define the G4PIIDATA environment variable as appropriate
  to get access to the PIXE data library (e.g. G4PII1.1)
  
  + To run without visualisation (batch mode):

    Go to $G4WORKDIR/bin/$G4SYSTEM
    Copy the file eRosita/eRosita.in to this directory. 
    The input file eRosita.in defines a simulation with 1000 protons 
    of energy 100 MeV.
    Start the simulation with: eRosita eRosita.in > eRosita.out

    The resulting files eRosita.out and the ASCII output file 
    TrackerPhotonEnergy.out are included in the example. Format
    and content of the output file are described below.

  + To run with visualisation:
    Go to $G4WORKDIR/bin/$G4SYSTEM
    Copy eRosita/vis.mac to to this directory.
    The macro file vis.mac calls the DAWN visualization driver to
    display the simulation of 100 protons with energy 100 MeV.

    An ASCII output file TrackerPhotonEnergy.out is created. However,
    this file may be empty in case the first 100 protons do not 
    produce any fluorescence photons that reach the tracker.
