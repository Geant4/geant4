// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TritonInelasticProcess.hh,v 1.3 2000-12-14 08:47:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Triton Inelastic Process
 // J.L. Chuma, TRIUMF, 25-Feb-1997
 // Last modified: 03-Apr-1997

// Class Description
// Process for Triton Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End
#ifndef G4TritonInelasticProcess_h
#define G4TritonInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4TritonInelasticProcess : public G4HadronicInelasticProcess
 class G4TritonInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4TritonInelasticProcess(
     const G4String& processName = "TritonInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4Triton::Triton() )
      G4HadronInelasticProcess( processName, G4Triton::Triton() )
    { }
        
    ~G4TritonInelasticProcess()
    { }
 };
 
#endif
 

