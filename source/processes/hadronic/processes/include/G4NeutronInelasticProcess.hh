// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronInelasticProcess.hh,v 1.3 2000-12-14 08:47:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Neutron Inelastic Process
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4NeutronInelasticProcess_h
#define G4NeutronInelasticProcess_h 1
 
// Class Description
// Process for Neutron Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4NeutronInelasticProcess : public G4HadronicInelasticProcess
 class G4NeutronInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4NeutronInelasticProcess(
     const G4String& processName = "NeutronInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4Neutron::Neutron() )
      G4HadronInelasticProcess( processName, G4Neutron::Neutron() )
    { }
    
    ~G4NeutronInelasticProcess()
    { }
 };
 
#endif
 
