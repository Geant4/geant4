// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPElastic.hh,v 1.2 1999-07-02 09:58:45 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: High Precision low E neutron tracking
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
 
#ifndef G4NeutronHPElastic_h
#define G4NeutronHPElastic_h 1

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"

class G4NeutronHPElastic : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPElastic();
  
  ~G4NeutronHPElastic();
  
  G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  G4int GetNiso() {return theElastic[0].GetNiso();}
  private:
  
  G4double * xSec;
  G4NeutronHPChannel * theElastic;
  G4String dirName;
  G4int numEle;
};

#endif
