// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPCaptureFS.hh,v 1.2 1999-06-29 18:43:49 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPCaptureFS_h
#define G4NeutronHPCaptureFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4ReactionProductVector.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPPhotonDist.hh"

class G4NeutronHPCaptureFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPCaptureFS();
  
  ~G4NeutronHPCaptureFS();
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPCaptureFS * theNew = new G4NeutronHPCaptureFS;
   return theNew;
  }
  
  private:
  
  G4double targetMass;
  
  G4NeutronHPPhotonDist theFinalStatePhotons;
  
  G4NeutronHPNames theNames;
};
#endif
