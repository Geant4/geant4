// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPLevel.hh"
#include "G4NeutronHPGamma.hh"

  G4NeutronHPLevel::G4NeutronHPLevel() 
  {
    nGammas = 0;
    theGammas = NULL;
  }

  G4NeutronHPLevel::~G4NeutronHPLevel() 
  {
    if(theGammas != NULL)
    {
      for(G4int i=0; i<nGammas; i++) delete theGammas[i];
    }
    delete [] theGammas;
  }

  void G4NeutronHPLevel::SetNumberOfGammas(G4int aGammas)
  {
    nGammas = aGammas;
    if(theGammas != NULL)
    {
      for(G4int i=0; i<nGammas; i++) delete theGammas[i];
    }
    delete [] theGammas; 
    theGammas = new G4NeutronHPGamma * [nGammas];
  }

  void G4NeutronHPLevel::SetGamma(G4int i, G4NeutronHPGamma * aGamma)
  {
    theGammas[i] = aGamma;
    SetLevelEnergy(aGamma->GetLevelEnergy());
  }

  G4double G4NeutronHPLevel::GetGammaEnergy(G4int i)
  {
    return theGammas[i]->GetGammaEnergy();
  }
  
  G4DynamicParticleVector * G4NeutronHPLevel::GetDecayGammas()
  {
    G4DynamicParticleVector * theResult;
    G4double sum = 0;
    G4double * running = new G4double[nGammas];
    running[0] = 0;
    G4int i;
    for(i=0; i<nGammas; i++)
    {
      if(i!=0) running[i]=running[i-1];
      running[i]+=theGammas[i]->GetWeight();
    }
    sum = running[nGammas-1];
    G4int it;
    G4double random = G4UniformRand();
    for(i=0; i<nGammas; i++)
    {
      it = i;
      if(random<running[i]/sum) break;
    }
    delete [] running;
    theResult = theGammas[it]->GetDecayGammas();
    return theResult;
  }
