// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AlphaInelasticProcess.hh,v 1.1 1999-01-07 16:13:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Alpha Inelastic Process
 // J.L. Chuma, TRIUMF, 21-Feb-1997
 // Last modified: 03-Apr-1997

#ifndef G4AlphaInelasticProcess_h
#define G4AlphaInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AlphaInelasticProcess : public G4HadronicInelasticProcess
 class G4AlphaInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AlphaInelasticProcess(
     const G4String& processName = "AlphaInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4Alpha::Alpha() )
      G4HadronInelasticProcess( processName, G4Alpha::Alpha() )
    { }
        
    ~G4AlphaInelasticProcess()
    { }
 };
 
#endif
 
