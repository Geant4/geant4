// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LambdaInelasticProcess.hh,v 1.1 1999-01-07 16:13:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process : Lambda Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997

 // Note:  there is no .cc file
 
#ifndef G4LambdaInelasticProcess_h
#define G4LambdaInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4LambdaInelasticProcess : public G4HadronicInelasticProcess
 class G4LambdaInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4LambdaInelasticProcess(
     const G4String& processName = "LambdaInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4Lambda::Lambda() )
      G4HadronInelasticProcess( processName, G4Lambda::Lambda() )
    { }
    
    ~G4LambdaInelasticProcess()
    { }
 };
 
#endif
 

