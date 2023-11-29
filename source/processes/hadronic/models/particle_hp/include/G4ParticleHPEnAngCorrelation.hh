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
#ifndef G4ParticleHPEnAngCorrelation_h
#define G4ParticleHPEnAngCorrelation_h 1

#include "G4Cache.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPProduct.hh"
#include "G4ParticleHPVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <fstream>

class G4ParticleDefinition;

class G4ParticleHPEnAngCorrelation
{
    struct toBeCached
    {
        G4ReactionProduct* theProjectileRP{nullptr};
        G4ReactionProduct* theTarget{nullptr};
        G4double theTotalMeanEnergy{-1.0};
        toBeCached() = default;
    };

  public:
    G4ParticleHPEnAngCorrelation()  // for G4ParticleHPCaptureFS::theMF6FinalState
    {
      theProjectile = G4Neutron::Neutron();
      theProducts = nullptr;
      inCharge = false;
      toBeCached val;
      fCache.Put(val);
      // theTotalMeanEnergy = -1.;
      fCache.Get().theTotalMeanEnergy = -1.;
      targetMass = 0.0;
      frameFlag = 0;
      nProducts = 0;
    }

    G4ParticleHPEnAngCorrelation(G4ParticleDefinition* proj) : theProjectile(proj)
    {
      theProducts = nullptr;
      inCharge = false;
      toBeCached val;
      fCache.Put(val);
      // theTotalMeanEnergy = -1.;
      fCache.Get().theTotalMeanEnergy = -1.;
      targetMass = 0.0;
      frameFlag = 0;
      nProducts = 0;
    }

    ~G4ParticleHPEnAngCorrelation() { delete[] theProducts; }

    inline void Init(std::istream& aDataFile)
    {
      inCharge = true;
      aDataFile >> targetMass >> frameFlag >> nProducts;
      theProducts = new G4ParticleHPProduct[nProducts];
      for (G4int i = 0; i < nProducts; i++) {
        theProducts[i].Init(aDataFile, theProjectile);
      }
    }

    G4ReactionProduct* SampleOne(G4double anEnergy);

    G4ReactionProductVector* Sample(G4double anEnergy);

    inline void SetTarget(G4ReactionProduct& aTarget)
    {
      fCache.Get().theTarget = &aTarget;
      for (G4int i = 0; i < nProducts; i++)
        theProducts[i].SetTarget(fCache.Get().theTarget);
    }

    inline void SetProjectileRP(G4ReactionProduct& aIncidentPart)
    {
      fCache.Get().theProjectileRP = &aIncidentPart;
      for (G4int i = 0; i < nProducts; i++)
        theProducts[i].SetProjectileRP(fCache.Get().theProjectileRP);
    }

    inline G4bool InCharge() { return inCharge; }

    inline G4double GetTargetMass() { return targetMass; }

    G4double GetTotalMeanEnergy()
    {
      return fCache.Get().theTotalMeanEnergy;  // cached in 'sample' call
    }

  private:
    // data members

    G4double targetMass;
    G4int frameFlag;  // =1: Target rest frame; =2: CMS system; incident in lab
    G4int nProducts;
    G4ParticleHPProduct* theProducts;
    G4bool inCharge;

    // cached values
    //
    G4Cache<toBeCached> fCache;

    G4ParticleDefinition* theProjectile;
};

#endif
