// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
 // Hadronic Process: Ion Inelastic Process
 // J.P. Wellisch, CERN, Apr. 14 2000
 // Last modified: 03-Apr-1997

#ifndef G4IonInelasticProcess_h
#define G4IonInelasticProcess_h 1
 
#include "G4HadronInelasticProcess.hh"
#include "G4GenericIon.hh"
 

 class G4IonInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4IonInelasticProcess(
     const G4String& processName = "IonInelastic" ) :

      G4HadronInelasticProcess( processName, G4GenericIon::GenericIon() )
    { }
        
    ~G4IonInelasticProcess()
    { }
 };

#endif

