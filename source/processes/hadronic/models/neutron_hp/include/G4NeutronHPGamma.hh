// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPGamma.hh,v 1.1 1999-01-07 16:13:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPGamma_h
#define G4NeutronHPGamma_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4DynamicParticleVector.hh"
#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4NeutronHPLevel.hh"

class G4NeutronHPGamma
{
  public:
  
  G4NeutronHPGamma() 
  {
    next = NULL;
  }
  ~G4NeutronHPGamma() {}
  
  G4bool Init(ifstream & aDataFile);
  
  inline void SetNext(G4NeutronHPLevel * aLevel)
  {
    next = aLevel;
  }
  
  G4DynamicParticleVector * GetDecayGammas()
  {
    G4DynamicParticleVector * theResult;
    if(next == NULL)
    {
      theResult = new G4DynamicParticleVector;
    }
    else
    {
      theResult = next->GetDecayGammas();
    }
    G4DynamicParticle * theNew = new G4DynamicParticle;
    theNew->SetDefinition(G4Gamma::Gamma());
    theNew->SetKineticEnergy(gammaEnergy);
    theResult->insert(theNew);
    return theResult;
  }
  
  inline G4double GetLevelEnergy()
  {
    return levelEnergy;
  }

  inline G4double GetGammaEnergy()
  {
    return gammaEnergy;
  }
  
  inline G4double GetWeight()
  {
    return probability;
  }

  private:
  
  G4double levelEnergy;
  G4double gammaEnergy;
  G4double probability;
  
  G4NeutronHPLevel * next;
};

#endif
