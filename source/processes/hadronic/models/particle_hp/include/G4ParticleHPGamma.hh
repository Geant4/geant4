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
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPGamma_h
#define G4ParticleHPGamma_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4DynamicParticleVector.hh"
#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4ParticleHPLevel.hh"

class G4ParticleHPGamma
{
  public:
  
  G4ParticleHPGamma();
  ~G4ParticleHPGamma();
  
  G4bool Init(std::istream & aDataFile);
  
  inline void SetNext(G4ParticleHPLevel * aLevel)
  {
    next = aLevel;
  }
  
  G4DynamicParticleVector * GetDecayGammas()
  {
    G4DynamicParticleVector * theResult;
    if(next == 0)
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
    theResult->push_back(theNew);
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
  
  G4ParticleHPLevel * next;
  static G4ThreadLocal int instancecount;
};

#endif
