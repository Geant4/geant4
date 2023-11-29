
        Geant4 extended examples - Hadronic processes
                      FlukaCern examples
        ----------------------------------------------

 Examples in this directory demonstrate how to make use of  
 the interface to `FLUKA` hadron inelastic physics in a G4 application.  
 
 The `FLUKA` interface itself is included,   
 and located in the `FlukaInterface` subdirectory.    
   
 There are 2 independent G4 examples at the process (interaction) level,   
 in order to study XS and final state respectively (in `ProcessLevel`),   
 and instructions at the physics list level (in `PhysicsListLevel`).  
 
 All examples benefit from `G4H1` extension (in `utils`),   
 allowing to make any G4 histogram also compatible with `Flair`.  
 
 Please directly follow `ProcessLevel/CrossSection/README.md`    
 and `ProcessLevel/FinalState/README.md`.  
 The `main` are included in `HadronNucleusXS.cc`   
 and `HadNucIneEvents.cc` respectively.    


ProcessLevel/CrossSection
------

This example allows the study of G4 cross-sections,   
and in addition, of the `FLUKA` hadron-nucleus reaction cross sections.  
The G4HadronicProcessStore is used to access cross sections.  
The resulting XS can be plotted via `ROOT` or `Flair`.  


ProcessLevel/FinalState/
------

This example shows how to simulate inelastic hadron-nucleus interactions,  
using G4 or `FLUKA` models.  
The resulting final states can be plotted via `ROOT` or `Flair`.  
Note that the Geant4 run-manager is not used.  


PhysicsListLevel
------

A `README.md` details how to select a physics list   
with hadron inelastic physics from FLUKA.
