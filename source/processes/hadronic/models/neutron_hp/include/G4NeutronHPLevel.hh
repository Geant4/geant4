// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPLevel.hh,v 1.1 1999-01-07 16:13:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPLevel_h
#define G4NeutronHPLevel_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4DynamicParticleVector.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
class G4NeutronHPGamma;

class G4NeutronHPLevel
{
  public:
  
  G4NeutronHPLevel() 
  {
    nGammas = 0;
    theGammas = NULL;
  }

  ~G4NeutronHPLevel();
  
  void SetNumberOfGammas(G4int aGammas);
  
  void SetGamma(G4int i, G4NeutronHPGamma * aGamma);
  
  G4DynamicParticleVector * GetDecayGammas();
    
  inline void SetLevelEnergy(G4double anEnergy)
  {
    levelEnergy = anEnergy;
  }
  
  inline G4double GetLevelEnergy()
  {
    return levelEnergy;
  }

  G4double GetGammaEnergy(G4int i);
  
  private:
  
  G4double levelEnergy;  

  G4int nGammas;
  G4NeutronHPGamma ** theGammas;
};

#endif
