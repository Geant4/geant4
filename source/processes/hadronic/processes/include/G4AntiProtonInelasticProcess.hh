// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiProtonInelasticProcess.hh,v 1.1 1999-01-07 16:13:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: AntiProton Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4AntiProtonInelasticProcess_h
#define G4AntiProtonInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AntiProtonInelasticProcess : public G4HadronicInelasticProcess
 class G4AntiProtonInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AntiProtonInelasticProcess(
     const G4String& processName = "AntiProtonInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4AntiProton::AntiProton() )
      G4HadronInelasticProcess( processName, G4AntiProton::AntiProton() )
    { }
    
    ~G4AntiProtonInelasticProcess()
    { }
 };
 
#endif
 

