// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiSigmaMinusInelasticProcess.hh,v 1.1 1999-01-07 16:13:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: AntiSigmaMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4AntiSigmaMinusInelasticProcess_h
#define G4AntiSigmaMinusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AntiSigmaMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4AntiSigmaMinusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AntiSigmaMinusInelasticProcess(
     const G4String& processName = "AntiSigmaMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4AntiSigmaMinus::AntiSigmaMinus() )
      G4HadronInelasticProcess( processName, G4AntiSigmaMinus::AntiSigmaMinus() )
    { }
    
    ~G4AntiSigmaMinusInelasticProcess()
    { }
 };
 
#endif
 

