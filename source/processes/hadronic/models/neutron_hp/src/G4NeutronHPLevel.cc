//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPLevel.hh"
#include "G4NeutronHPGamma.hh"

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
    G4int it(0);
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
