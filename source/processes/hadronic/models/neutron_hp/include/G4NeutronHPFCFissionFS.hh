// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFCFissionFS.hh,v 1.1 1999-01-07 16:12:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFCFissionFS_h
#define G4NeutronHPFCFissionFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFissionBaseFS.hh"

class G4NeutronHPFCFissionFS : public G4NeutronHPFissionBaseFS
{
  public:
  
  G4NeutronHPFCFissionFS(){ hasXsec = false; }
  ~G4NeutronHPFCFissionFS(){}
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4DynamicParticleVector * ApplyYourself(G4int nNeutrons);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPFCFissionFS * theNew = new G4NeutronHPFCFissionFS;
   return theNew;
  }
  
  private:
  G4ParticleChange * ApplyYourself(const G4Track & theTrack) { return NULL; }
    
};
#endif
