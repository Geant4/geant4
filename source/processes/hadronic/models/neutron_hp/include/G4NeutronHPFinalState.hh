// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFinalState.hh,v 1.2 1999-06-29 18:43:59 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFinalState_h
#define G4NeutronHPFinalState_h

#include "G4Material.hh"
#include "G4FastVector.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPVector.hh"

class G4NeutronHPFinalState
{
public:

  G4NeutronHPFinalState();
  
  virtual ~G4NeutronHPFinalState();

  virtual void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType) = 0;
  virtual G4ParticleChange * ApplyYourself(const G4Track & theTrack) 
  {
    G4Exception("G4ParticleChange * ApplyYourself(const G4Track & theTrack) needs implementation");
    return NULL;
  }
  
  // of course this would better be Done templating G4NeutronHPChannel..., 
  // but due to peculiarities of the DEC compiler, this way it
  // is easier to maintain.
  virtual G4NeutronHPFinalState * New() = 0;
  
  G4bool HasXsec() {return hasXsec;}
  G4bool HasFSData() {return hasFSData;}
  G4bool HasAnyData() {return hasAnyData;}
  
  virtual G4double GetXsec(G4double anEnergy) { return 0; }
  virtual G4NeutronHPVector * GetXsec() { return NULL; }
  
  G4double GetZ() { return theBaseZ; }
  G4double GetN() { return theBaseA; }
  
  protected:
  
  G4bool hasXsec;
  G4bool hasFSData;
  G4bool hasAnyData;
  G4NeutronHPNames theNames;
  
  G4ParticleChange theResult;
  
  G4double theBaseA;
  G4double theBaseZ;

  private:
  
};
#endif
