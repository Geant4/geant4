// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KaonMinusInelasticProcess.hh,v 1.2 1999-12-15 14:53:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: KaonMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 12-Feb-1997
 // Last modified: 03-Apr-1997

#ifndef G4KaonMinusInelasticProcess_h
#define G4KaonMinusInelasticProcess_h 1
 
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
 
