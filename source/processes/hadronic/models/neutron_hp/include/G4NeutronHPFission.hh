// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFission.hh,v 1.2 1999-07-02 09:59:01 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: High Precision low E neutron tracking
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
  
#ifndef G4NeutronHPFission_h
#define G4NeutronHPFission_h 1

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"

#include "G4NeutronHPFissionFS.hh"

class G4NeutronHPFission : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPFission();

  ~G4NeutronHPFission();
  
  G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
  
  G4NeutronHPFissionFS theFS;
  
  private:
  
  G4double * xSec;
  G4NeutronHPChannel * theFission;
  G4String dirName;
  G4int numEle;
  static G4String theNames[3];
};

#endif
