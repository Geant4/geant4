// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiXiZeroInelasticProcess.hh,v 1.3 2000-12-14 08:47:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: AntiXiZero Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4AntiXiZeroInelasticProcess_h
#define G4AntiXiZeroInelasticProcess_h 1
 
// Class Description
// Process for AntiXiZero Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AntiXiZeroInelasticProcess : public G4HadronicInelasticProcess
 class G4AntiXiZeroInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AntiXiZeroInelasticProcess(
     const G4String& processName = "AntiXiZeroInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4AntiXiZero::AntiXiZero() )
      G4HadronInelasticProcess( processName, G4AntiXiZero::AntiXiZero() )
    { }
    
    ~G4AntiXiZeroInelasticProcess()
    { }
 };
 
#endif
 

