// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KaonPlusInelasticProcess.hh,v 1.2 1999-12-15 14:53:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: KaonPlus Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Feb-1997
 // Last modified: 03-Apr-1997

#ifndef G4KaonPlusInelasticProcess_h
#define G4KaonPlusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4KaonPlusInelasticProcess : public G4HadronicInelasticProcess
 class G4KaonPlusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4KaonPlusInelasticProcess(
     const G4String& processName = "KaonPlusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4KaonPlus::KaonPlus() )
      G4HadronInelasticProcess( processName, G4KaonPlus::KaonPlus() )
    { }
    
    ~G4KaonPlusInelasticProcess()
    { }
 };
 
#endif
 
