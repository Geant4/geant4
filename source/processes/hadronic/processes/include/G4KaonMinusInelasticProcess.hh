// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KaonMinusInelasticProcess.hh,v 1.3 2000-12-14 08:47:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: KaonMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 12-Feb-1997
 // Last modified: 03-Apr-1997

#ifndef G4KaonMinusInelasticProcess_h
#define G4KaonMinusInelasticProcess_h 1
 
// Class Description
// Process for KaonMinus Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4KaonMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4KaonMinusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4KaonMinusInelasticProcess(
     const G4String& processName = "KaonMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4KaonMinus::KaonMinus() )
      G4HadronInelasticProcess( processName, G4KaonMinus::KaonMinus() )
    { }
    
    ~G4KaonMinusInelasticProcess()
    { }
 };
 
#endif
 
