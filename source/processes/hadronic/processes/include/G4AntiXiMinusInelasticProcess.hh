// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiXiMinusInelasticProcess.hh,v 1.1 1999-01-07 16:13:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: AntiXiMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4AntiXiMinusInelasticProcess_h
#define G4AntiXiMinusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AntiXiMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4AntiXiMinusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AntiXiMinusInelasticProcess(
     const G4String& processName = "AntiXiMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4AntiXiMinus::AntiXiMinus() )
      G4HadronInelasticProcess( processName, G4AntiXiMinus::AntiXiMinus() )
    { }
    
    ~G4AntiXiMinusInelasticProcess()
    { }
    
 };
 
#endif
 

