// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4XiMinusInelasticProcess.hh,v 1.3 2000-12-14 08:47:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // G4 Process: XiMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997

// Class Description
// Process for XiMinus Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

#ifndef G4XiMinusInelasticProcess_h
#define G4XiMinusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4XiMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4XiMinusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4XiMinusInelasticProcess(
     const G4String& processName = "XiMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4XiMinus::XiMinus() )
      G4HadronInelasticProcess( processName, G4XiMinus::XiMinus() )
    { }
    
    ~G4XiMinusInelasticProcess()
    { }
 };
 
#endif
