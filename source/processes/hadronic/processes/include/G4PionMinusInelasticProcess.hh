// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PionMinusInelasticProcess.hh,v 1.2 1999-12-15 14:53:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // PionMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997

#ifndef G4PionMinusInelasticProcess_h
#define G4PionMinusInelasticProcess_h 1
 
//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4PionMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4PionMinusInelasticProcess : public G4HadronInelasticProcess
 {
    
 public:
    
    G4PionMinusInelasticProcess(
     const G4String& processName = "PionMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4PionMinus::PionMinus() )
      G4HadronInelasticProcess( processName, G4PionMinus::PionMinus() )
    { }
    
    ~G4PionMinusInelasticProcess()
    { }
    
 };
 
#endif
