// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PionPlusInelasticProcess.hh,v 1.1 1999-01-07 16:13:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // PionPlus Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997
 
#ifndef G4PionPlusInelasticProcess_h
#define G4PionPlusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4PionPlusInelasticProcess : public G4HadronicInelasticProcess
 class G4PionPlusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4PionPlusInelasticProcess(
     const G4String& processName = "PionPlusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4PionPlus::PionPlus() )
      G4HadronInelasticProcess( processName, G4PionPlus::PionPlus() )
    { }
    
    ~G4PionPlusInelasticProcess()
    { }
    
 };
 
#endif
