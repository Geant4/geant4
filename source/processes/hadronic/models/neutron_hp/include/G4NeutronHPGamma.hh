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
//
// $Id: G4NeutronHPGamma.hh,v 1.12 2005/06/04 13:44:43 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#ifndef G4NeutronHPGamma_h
#define G4NeutronHPGamma_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4DynamicParticleVector.hh"
#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4NeutronHPLevel.hh"

class G4NeutronHPGamma
{
  public:
  
  G4NeutronHPGamma();
  ~G4NeutronHPGamma();
  
  G4bool Init(std::ifstream & aDataFile);
  
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
  
  G4NeutronHPLevel * next;
  static int instancecount;
};

#endif
