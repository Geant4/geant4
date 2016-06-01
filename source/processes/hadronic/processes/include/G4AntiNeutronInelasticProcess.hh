// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiNeutronInelasticProcess.hh,v 2.0 1998/07/02 16:36:09 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
 // Hadronic Process: AntiNeutron Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4AntiNeutronInelasticProcess_h
#define G4AntiNeutronInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AntiNeutronInelasticProcess : public G4HadronicInelasticProcess
 class G4AntiNeutronInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AntiNeutronInelasticProcess(
     const G4String& processName = "AntiNeutronInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4AntiNeutron::AntiNeutron() )
      G4HadronInelasticProcess( processName, G4AntiNeutron::AntiNeutron() )
    { }
    
    ~G4AntiNeutronInelasticProcess()
    { }
 };
 
#endif
 

