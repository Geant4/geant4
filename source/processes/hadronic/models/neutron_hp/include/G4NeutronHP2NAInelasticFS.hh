// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHP2NAInelasticFS.hh,v 1.2 1999-06-29 18:43:44 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHP2NAInelasticFS_h
#define G4NeutronHP2NAInelasticFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPInelasticBaseFS.hh"
#include "G4NeutronHPAngular.hh"
#include "G4NeutronHPEnergyDistribution.hh"
#include "G4NeutronHPEnAngCorrelation.hh"
#include "G4NeutronHPPhotonDist.hh"

class G4NeutronHP2NAInelasticFS : public G4NeutronHPInelasticBaseFS
{
  public:
  
  G4NeutronHP2NAInelasticFS();
  ~G4NeutronHP2NAInelasticFS();
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHP2NAInelasticFS * theNew = new G4NeutronHP2NAInelasticFS;
   return theNew;
  }
  
  private:
  
};
#endif
