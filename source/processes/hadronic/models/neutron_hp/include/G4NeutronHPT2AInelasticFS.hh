// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPT2AInelasticFS.hh,v 1.2 1999-06-29 18:44:15 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPT2AInelasticFS_h
#define G4NeutronHPT2AInelasticFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPInelasticBaseFS.hh"

class G4NeutronHPT2AInelasticFS : public G4NeutronHPInelasticBaseFS
{
  public:
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPT2AInelasticFS * theNew = new G4NeutronHPT2AInelasticFS;
   return theNew;
  }
  
  private:
};
#endif
