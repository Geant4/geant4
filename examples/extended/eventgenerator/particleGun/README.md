\page ExampleparticleGun Example particleGun 
			    
 History:
 - 10-06-2010 : Makoto Asai - merge into one example
 - 13-05-2010 : Michel Maire - create as three examples
  
 This example demonstrates 5 ways of the use of G4ParticleGun shooting
 primary particles in different cases. These are
 -#  uniform particle direction in a given solid angle
 -#  Generate several vertices and particles per event
 -#  Show how to sample a tabulated function (eg. energy spectrum)  
 -#  Divergent beam in an arbitrary direction 
 -#  Shooting primaries in spherical coordinates with rotation matrix.

 These use cases can be chosen by a UI command
```
/gunExample/selectGunAction actionID
```
 where <i>actionID</i> corresponds to the above cases.
 	
## Geometry construction

 It is a simple box which represents an 'infinite' homogeneous medium.
  
## Physics list

PhysicsList defines only geantino and transportation process.
         	
## Primary generator

 There are 5 concrete primary generator action classes
 (PrimaryGeneratorActionN, N=0,1,2,3,4) which can be used as an independent sample
 code. 
 PrimaryGeneratorAction is the class which uses and switches between these 
 5 concrete action classes. Each concrete generator action shoots geantinoes
 in a distribution decribed below.
 
### 0. Uniform particle direction in a given solid angle

 spherical angles (alpha,psi) respective to z axis
 Histograms 5,6 show momentum direction in master frame.

 See: PrimaryGeneratorAction0::GeneratePrimaries()
 
### 1. Generate several vertices and particles per event

 - particle 1 : a geantino uniformly randomized on a cylinder surface.
 - particle 2 and 3 : symetric to particle 1.
 In addition, time_zero of each event is randomized.
 
 See: PrimaryGeneratorAction1::GeneratePrimaries()

### 2. Show how to sample a tabulated function (energy spectrum)

 Energy is sampled from a tabulated function defined in InitFunction().
 The function is assumed positive, linear per segment, continuous.
 Two sampling methods are illustrated : RejectAccept() and InverseCumul()
 (see Particle Data Group: pdg.lbl.gov --> Monte Carlo techniques)

 Histogram 1 shows generated energy spectrum.

 See: PrimaryGeneratorAction2::GeneratePrimaries()
       
### 3. Divergent beam in an arbitrary direction with rotation matrix

 A geantino uniformly randomized around a given direction (theta, phi). 
 One wants to limit particle direction uniformly around this direction.
 First, one generates momentum direction in the master frame (eg. World).
 AlphaMax = opening angle around z axis.
 Then one rotates momentum in local frame, using rotateUz() function.
 RotateUz() transforms uz to newUz. It is composition of two simple rotations:
 theta around oy, then phi around oz (non commutative). \n
 See: https://proj-clhep.web.cern.ch/proj-clhep/  \n
 Histograms 5,6 show momentum direction in local frame.
  
 See: PrimaryGeneratorAction3::GeneratePrimaries()

### 4. Shooting primaries in spherical coordinates with rotation matrix

 A geantino uniformly randomized within a spherical shell.

 a) Vertex position
 One wishes to shoot uniformly within a spherical shell.
 One works in spherical coordinates. One uses inverse cumulative method with
 analytical formulae. \n
 Histograms 2,3,4 demonstrate uniform distribution of vertex position.

 b) Momentum direction
 One wants to limit particle direction uniformly within (alphaMin, alphaMax).
 First, one generates momentum direction in the master frame (eg. World).
 Then, one rotates momentum in vertex_position frame, using rotateUz() function.
 RotateUz() transforms uz to ur. It is composition of two elementary rotations:
 theta around oy, then phi around oz (non commutative). \n
 See: https://proj-clhep.web.cern.ch/proj-clhep/ \n
 
 Histograms 5,6 show momentum direction in vertex_position frame.

 See: PrimaryGeneratorAction4::GeneratePrimaries()

## Visualisation
 
 Visualization Manager is set in the main().
 Initialisation of the drawing is done via the commands
 `/vis/..`` in the macro vis.mac. This macro is automatically read from the main 
 in case of interactive running mode.
 
## How to start ?
 		
  - Execute particleGun in 'batch' mode from macro files
```
% ./particleGun   run1.mac
```
 		
  - Execute particleGun2 in 'interactive mode' with visualization
```
% ./particleGun
....
Idle>      ---> type your commands. For instance:
Idle> /gunExample/selectGunAction 1
Idle> /run/beamOn 10
....
Idle> exit
```

  - Macros provided in this example:
      - run0.mac: see above primaryGenerator, case 0
      - run1.mac: see above primaryGenerator, case 1
      - run2.mac: see above primaryGenerator, case 2
      - run3.mac: see above primaryGenerator, case 3
      - run4.mac: see above primaryGenerator, case 4
      - particleGun.in: macro used in Geant4 testing (similar to run0.mac)
	
## Histograms
 
  particleGun produces several 1D histograms which are saved as
  particleGun.root by default.

    - 1 : energy spectrum dN/dE = f(E)
    - 2 : vertex position: radial distr dN/dv = f(r)
    - 3 : vertex position: cos(theta)
    - 4 : vertex position: phi
    - 5 : particle direction in local frame: cos(alpha)
    - 6 : particle direction in local frame: psi      

  Please note that histogram 1 will be filled only if you use 
  PrimaryGeneratorAction2, histos 5,6 with PrimaryGeneratorAction0 and 3
  and 2 through 6 will be filled with PrimaryGeneratorAction4.

  The histograms are managed by the HistoManager class and its messenger. 
  The histos can be individually activated with the command :
```
/analysis/h1/set id nbBins  valMin valMax unit 
```
  where unit is the desired unit for the histo (MeV or keV, deg or mrad, etc..)
   
  One can control the name of the histograms file with the command:
```
/analysis/h1/setFileName  name  (default particleGun)
```
   
   It is possible to choose the format of the histogram file : root (default),
   xml, csv by setting default file type in HistoManager::Book().
   
   It is also possible to print selected histograms on an ascii file:
```
/analysis/h1/setAscii id
```
   All selected histos will be written on a file name.ascii (default gunExample)
