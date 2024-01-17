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
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//
#ifndef G4ParticleHPProduct_h
#define G4ParticleHPProduct_h 1

#include "G4Cache.hh"
#include "G4ReactionProductVector.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ParticleHPVector.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>

class G4ParticleDefinition;

enum G4HPMultiMethod
{
  G4HPMultiPoisson,
  G4HPMultiBetweenInts
};

class G4ParticleHPProduct
{
    struct toBeCached
    {
        G4ReactionProduct* theProjectileRP{nullptr};
        G4ReactionProduct* theTarget{nullptr};
        G4int theCurrentMultiplicity{-1};
        toBeCached() = default;
    };

  public:
    G4ParticleHPProduct();
    ~G4ParticleHPProduct();

    void Init(std::istream& aDataFile, const G4ParticleDefinition* projectile);

    G4int GetMultiplicity(G4double anEnergy);
    G4ReactionProductVector* Sample(G4double anEnergy, G4int nParticles);

    G4double GetMeanYield(G4double anEnergy) { return theYield.GetY(anEnergy); }

    void SetProjectileRP(G4ReactionProduct* aIncidentPart)
    {
      fCache.Get().theProjectileRP = aIncidentPart;
    }

    void SetTarget(G4ReactionProduct* aTarget) { fCache.Get().theTarget = aTarget; }

    inline G4ReactionProduct* GetTarget() { return fCache.Get().theTarget; }

    inline G4ReactionProduct* GetProjectileRP() { return fCache.Get().theProjectileRP; }

    inline G4double MeanEnergyOfThisInteraction()
    {
      G4double result = 0.0;
      if (theDist != nullptr) {
        result = theDist->MeanEnergyOfThisInteraction();
        result *= fCache.Get().theCurrentMultiplicity;
      }
      return result;
    }

    inline G4double GetQValue() { return theActualStateQValue; }

    // TK120515 For migration of frameFlag (MF6 LCT) = 3 in
    // G4ParticleHPEnAngCorrelation
    G4double GetMassCode() { return theMassCode; }
    G4double GetMass() { return theMass; }

  private:
    G4double theMassCode{0.0};
    G4double theMass{0.0};
    G4double theGroundStateQValue{0.0};
    G4double theActualStateQValue{0.0};
    G4int theIsomerFlag{0};
    G4int theDistLaw{-1};  // redundant
    G4VParticleHPEnergyAngular* theDist{nullptr};

    // cashed values
    //
    G4Cache<toBeCached> fCache;

    G4HPMultiMethod theMultiplicityMethod;
    G4ParticleHPVector theYield;
};

#endif
