// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPCapture.hh,v 1.3 2000-08-10 12:52:58 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
 
#ifndef G4NeutronHPCapture_h
#define G4NeutronHPCapture_h 1

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"

class G4NeutronHPCapture : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPCapture();

  ~G4NeutronHPCapture();

  G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
  
  G4double * xSec;
  G4NeutronHPChannel * theCapture;
  G4String dirName;
  G4int numEle;
  G4int it;
  
  G4ParticleChange theResult;
};

#endif
