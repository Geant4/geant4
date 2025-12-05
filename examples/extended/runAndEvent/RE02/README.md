\page ExampleRE02 Example RE02

 This example simulates a simplified water phantom measurement
 in medical application with demonstration of primitive scorers.
 This example also demonstrates nested parameterised volume which
 realizes segmented boxes using a combination of replicated volumes
 and a parameterised volume.

  ---- (Tips)

 This example creates 100 x 100 x 200 boxes using Nested Parameterised
 Volume for realistic situation of medical application. 
 This is very memory consumption if normal Parameterised Volume is used, 
 and needs roughly more than 1 GB memory for execution. However, 
 NestedParameterised volume effectively works to reduce the memory consumption,
 and it only needs less than 100 MB memory for execution.

## GEOMETRY DEFINITION
 
  The setup contains a water phantom as target by default. The world volume
  is 200 cm x 200 cm x 200 cm box filled with air. The water phantom is box shape 
  and the size of 200 mm x 200 mm x 400 mm. The volume of water phantom is divided
  into 100 x 100 x 1 towers using replicated volume,(RE02DetectorConstruction), 
  and then those towers are segmented into 200 boxes with respect to z axis 
  using nested parameterized volume,(RE02NestedPhantomParameterisation).
  e.g. The volume of water phantom is divided into 100 x 100 x 200 boxes,
  and a voxel size is 2.0 mm x 2.0 mm x 2.0 mm.

  For demonstration purpose of the nested parameterised volume,
  (RE02NestedPhantomParameterisation), materials are assigned as water (lead) 
  in even (odd) order segments, alternately. 
  The simulation for homogeneous water phantom is also possible using an option.

  ---- Tips(1)
  
  If you want to reduce number of segments of water phantom,
  please change following numbers which represent number of segments
  in x, y, z axis, respectively.The following code can be found in
  RE02.cc.
```cpp
  RE02DetectorConstruction* detector = new RE02DetectorConstruction;
  detector->SetNumberOfSegmentsInPhantom(100, 100, 200);
  //                                      Nx,  Ny,  Nz
```
  ---- Tips(2)
  
  If you want to set all materials to water,
   please use the following method. The following code can be found in
  RE02.cc.
```cpp
   detector->SetLeadSegment(false); // Homogeneous water phantom

```

  The geometry and sensitive detector are constructed in 
  RE02DetectorConstruction class.
  (See the SCORER section for detail descriptions about sensitive detector.)
         
## PHYSICS LIST
 
  This example uses the QGS_BIC physics list.
 	
## RUNS and EVENTS
 
### Primary particles
     
  The primary kinematics consists of a single particle which hits the
  target perpendicular to the input face. The default type of the particle
  and its energy are set in the RE02PrimaryGeneratorAction class.
  However it can be changed via the G4 build-in commands of ParticleGun 
  class. 
   The RE02PrimaryGeneratorAction class introduces a beam spot size 
  that makes initial particle position of x,y randomized using a Gaussian
  random function, where the center position is fixed to (0,0). 
  The standard deviation of the beam spot size is given in 
  RE02PrimaryGeneratorAction as 10 mm.

### Event
  
  An EVENT represents a simulation of one primary particle.
  A RUN is a set of events.
  
  The user has control:
     - at Begin and End of each run   (class RunAction)
     - at Begin and End of each event (class EventAction)
     - at Begin and End of each track (class TrackingAction, not used here)
     - at End of each step (class SteppingAction, not used here)
	
## SCORER

### Concrete Scorer

   This example introduces concrete primitive scorer (PS) and filter 
   classes for easy scoring. Those primitive scorers are registered to
   G4MultiFunctionalDetector which is a concrete class of sensitive 
   detector(SD). Then the G4MultiFunctionalDetector is attached to 
   the logical volume of sensitive geometry.
   A MultiFunctionalDetector, PrimitiveScorers, and SDFilters are
   created and assigned to the logical volume of water phantom in 
   DetectorConstruction.

   A primitive scorer can score one kind of physical quantity, and
   creates one hits collection per event. The quantity is collected in
   G4THitsMap with the copy number of geometry. Here collection name is 
   given as "MultiFunctionalDetector Name"/"PrimitiveScorer Name".
   A primitive scorer can have one filter (SDFilter) for selecting hits 
   to be used for the quantity.

   Since the geometry is constructed using nested parameterisation,
   the copy number of geometry is defined as follows,
```
  copy number of geometry =  iy*Nx*Ny+ix*Nz+iz,
```

   where Nx,Ny,Nz is total number of segmentation in x, y, and z axis,respectively,
   and ix,iy,iz is a copy number of the mother volume, the grand mother volume, 
   and this volume, respectively.
    This conversion is described in GetIndex() method in PrimitiveScorer.

### The physical quantities scored in this example are:

   - Total energy deposit \n
      - unit: Energy,                  collName: totalEDep
   - Energy deposit by protons \n
      - unit: Energy,                  collName: protonEDep
   - Number of steps of protons \n
      - unit: - ,                      collName: protonNStep   
   - Cell Flux of charged tracks which pass through the geometry\n
      - unit: Length/Volume,          collName:  chargedPassCellFlux
   - Cell Flux of all charged tracks\n
      - unit: Length/Volume,           collName:  chargedCellFlux      
   - Flux of charged particle at -Z surface of the BOX geometry,
     where incident angle at the surface is taken into account.\n
      - unit: Surface^(-1),            collName:  chargedSurfFlux 
   - Surface current of gamma at -Z surface of the BOX geometry. 
     The energy of gammas are from    1. keV  to   10. keV.
     The incident angle is not taken into account.\n
      - unit: Surface^(-1),            collName:  gammaSurfCurr000
   - Same as previous one, but different energy bin.
     The energy of gammas are from   10. keV  to  100. keV.\n
      - unit: Surface^(-1),           collName:  gammaSurfCurr001
   - Same as previous one, but different energy bin.
     The energy of gammas are from  100. keV  to    1. MeV. \n
      - unit: Surface^(-1),           collName:  gammaSurfCurr002
   - Same as previous one, except for energy bin.   
     The energy of gammas are from    1. MeV  to   10. MeV. \n
      - unit: Surface^(-1),            collName:  gammaSurfCurr003
 
### Accumulating quantities during a RUN
     
  A PrimitiveScorer creates one hits collection per event.
  The physical quantity in the hits collection need to be accumulated
  into another G4THitsMap object during a RUN, in order to obtain
  integrated flux or dose in a RUN. The accumulation of quantities 
  are done at RE02Run class.
 
  RE02Run class can automatically generate G4THitsMap objects for a RUN,
  and accumulate physical quantities of an event into it. The accumulation
  is done at RE02Run::RecordEvent().
  
### Generate a Run object, and print results
  
   The RE02Run object is generated at RE02RunAction::GenerateRun().      
  The accumulated physical quantities are printed at the end of RUN
  ( RE02RunAction::EndOfRunAction() ). This example prints only selected
  physical quantities.

 				
## VISUALIZATION
 
  The Visualization Manager is set in the main().
  The initialization of the drawing is done via a set of /vis/ commands
  in the macro vis.mac. This macro is automatically read from 
  the main when running in interactive mode.
 
  The tracks are automatically drawn at the end of event and erased at 
  the beginning of the next run.

  The visualization (with OpenGL driver) assumes two things:
     -# the visualization & interfaces categories have been compiled
            with the environment variable G4VIS_BUILD_OPENGLX_DRIVER.
     -# RE02.cc has been compiled with G4VIS_USE_OPENGLX.   

  (The same with DAWNFILE instead of OPENGLX)
     
     
## USER INTERFACES
  
  The default command interface, called G4UIterminal, is done via
  standard G4cin/G4cout.
  On Linux and Sun-cc on can use a smarter command interface G4UItcsh.
  It is enough to set the environment variable G4UI_USE_TCSH before
  compiling RE02.cc
 
       	
## HOW TO START ?
 
  - Execute RE02 in 'batch' mode from macro files (without visualization)
```
% ./RE02   run1.mac
```
 		
  - Execute RE02 in 'interactive mode' with visualization
```
% ./RE02
....
Idle> type your commands. For instance:
Idle> /run/beamOn 10
....
Idle> /control/execute run2.mac
....
Idle> exit
```
     
  - Macros are for different primary particles.
    - vis.mac  :  200 MeV proton with visualization
    - run1.mac :  150 MeV proton
    - run2.mac :  195 MeV/u Carbon ion
    - run3.mac :   30 MeV electron
    - run4.mac :   60 keV gamma

 	
      
