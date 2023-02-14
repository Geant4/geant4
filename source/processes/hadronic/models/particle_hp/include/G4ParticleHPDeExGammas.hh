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
#ifndef G4ParticleHPDeExGammas_h
#define G4ParticleHPDeExGammas_h 1

#include <fstream>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"
#include "G4ReactionProductVector.hh"
#include "G4Gamma.hh"
#include "G4ParticleHPLevel.hh"
#include "G4ParticleHPGamma.hh"
#include "G4ReactionProduct.hh"

class G4ParticleHPDeExGammas
{
  public:
  
  G4ParticleHPDeExGammas() 
  {
  }
  ~G4ParticleHPDeExGammas() 
  {
    delete [] levelStart;
    delete [] levelSize;
    delete [] theLevels;
  }
  
  void Init(std::istream & aDataFile);

  inline G4ReactionProductVector * GetDecayGammas(G4int aLevel)
  {
    if(aLevel>nLevels-1 || aLevel<0) return nullptr;
    if(nLevels==0) return new G4ReactionProductVector();
    G4ReactionProductVector * result = new G4ReactionProductVector;
    G4DynamicParticleVector * theResult;

    theResult = theLevels[aLevel]. GetDecayGammas();
    G4ReactionProduct * theCurrent;
    for(unsigned int i=0; i<theResult->size(); ++i)
    {
      theCurrent = new G4ReactionProduct;
      *theCurrent = *(theResult->operator[](i));
      delete theResult->operator[](i);
      G4double costheta = 2.*G4UniformRand()-1;
      G4double theta = std::acos(costheta);
      G4double phi = CLHEP::twopi*G4UniformRand();
      G4double sinth = std::sin(theta);
      G4double en = theCurrent->GetTotalMomentum();
      G4ThreeVector temp(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*costheta );
      theCurrent->SetMomentum( temp ) ;
      result->push_back(theCurrent);
    }
    delete theResult;
    return result;
  }
  
  inline G4ParticleHPLevel * GetLevel(G4int i)
  {
    if(i>nLevels-1) return nullptr;
    return theLevels+i;
  }
  
  inline G4int GetNumberOfLevels() { return nLevels; }
  
  inline G4double GetLevelEnergy(G4int aLevel)
  {
    if(aLevel>nLevels-1 || aLevel<0) return 0;
    G4double result = theLevels[aLevel].GetLevelEnergy();
    return result;
  }
  private:
  
  G4int * levelStart = nullptr;
  G4int * levelSize = nullptr;
  G4int nLevels = 0;
  G4ParticleHPLevel * theLevels = nullptr;
};

#endif
