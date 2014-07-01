//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4ParticleHPLevel.hh"
#include "G4ParticleHPGamma.hh"

  G4ParticleHPLevel::~G4ParticleHPLevel() 
  {
    if(theGammas != 0)
    {
      for(G4int i=0; i<nGammas; i++) delete theGammas[i];
    }
    delete [] theGammas;
  }

  void G4ParticleHPLevel::SetNumberOfGammas(G4int aGammas)
  {
    nGammas = aGammas;
    if(theGammas != 0)
    {
      for(G4int i=0; i<nGammas; i++) delete theGammas[i];
    }
    delete [] theGammas; 
    theGammas = new G4ParticleHPGamma * [nGammas];
  }

  void G4ParticleHPLevel::SetGamma(G4int i, G4ParticleHPGamma * aGamma)
  {
    theGammas[i] = aGamma;
    SetLevelEnergy(aGamma->GetLevelEnergy());
  }

  G4double G4ParticleHPLevel::GetGammaEnergy(G4int i)
  {
    return theGammas[i]->GetGammaEnergy();
  }
  
  G4DynamicParticleVector * G4ParticleHPLevel::GetDecayGammas()
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
      if(random*sum < running[i]) break;
    }
    delete [] running;
    theResult = theGammas[it]->GetDecayGammas();
    return theResult;
  }
