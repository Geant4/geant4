// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFissionFS.hh,v 1.2 1999-06-29 18:44:00 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFissionFS_h
#define G4NeutronHPFissionFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"

#include "G4NeutronHPFCFissionFS.hh"
#include "G4NeutronHPSCFissionFS.hh"
#include "G4NeutronHPTCFissionFS.hh"
#include "G4NeutronHPLCFissionFS.hh"
#include "G4NeutronHPFSFissionFS.hh"

class G4NeutronHPFissionFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPFissionFS();
  ~G4NeutronHPFissionFS();
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPFissionFS * theNew = new G4NeutronHPFissionFS;
   return theNew;
  }
        
  private:
  
  G4NeutronHPFSFissionFS theFS;
  G4NeutronHPFCFissionFS theFC;
  G4NeutronHPSCFissionFS theSC;
  G4NeutronHPTCFissionFS theTC;
  G4NeutronHPLCFissionFS theLC;
    
};
#endif
