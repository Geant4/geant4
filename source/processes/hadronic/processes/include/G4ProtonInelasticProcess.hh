// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProtonInelasticProcess.hh,v 1.2 1999-12-15 14:53:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Proton Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997
 //
 // Note:  there is no .cc file
 
#ifndef G4ProtonInelasticProcess_h
#define G4ProtonInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4ProtonInelasticProcess : public G4HadronicInelasticProcess
 class G4ProtonInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4ProtonInelasticProcess(
     const G4String& processName = "ProtonInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4Proton::Proton() )
      G4HadronInelasticProcess( processName, G4Proton::Proton() )
    { }
    
    ~G4ProtonInelasticProcess()
    { }
    
 };
 
#endif
 
