// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SigmaMinusInelasticProcess.hh,v 1.2 1999-12-15 14:53:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // G4 Process: SigmaMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997

#ifndef G4SigmaMinusInelasticProcess_h
#define G4SigmaMinusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4SigmaMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4SigmaMinusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4SigmaMinusInelasticProcess(
     const G4String& processName = "SigmaMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4SigmaMinus::SigmaMinus() )
      G4HadronInelasticProcess( processName, G4SigmaMinus::SigmaMinus() )
    { }
    
    ~G4SigmaMinusInelasticProcess()
    { }
 };
 
#endif
