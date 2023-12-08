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
//
// 080721 Add ClearHistories() method by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//
#ifndef G4ParticleHPContEnergyAngular_h
#define G4ParticleHPContEnergyAngular_h 1

#include "G4Cache.hh"
#include "G4InterpolationManager.hh"
#include "G4ParticleHPContAngularPar.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>

class G4ParticleDefinition;

class G4ParticleHPContEnergyAngular : public G4VParticleHPEnergyAngular
{
public:

  G4ParticleHPContEnergyAngular(const G4ParticleDefinition* proj);

  ~G4ParticleHPContEnergyAngular() override;

  void Init(std::istream& aDataFile) override;

  G4double MeanEnergyOfThisInteraction() override;
  G4ReactionProduct* Sample(G4double anEnergy, G4double massCode, G4double mass) override;
  void ClearHistories() override;

  G4ParticleHPContEnergyAngular(G4ParticleHPContEnergyAngular&) = delete;
  G4ParticleHPContEnergyAngular& operator=
  (const G4ParticleHPContEnergyAngular& right) = delete;

private:
  G4double theTargetCode{-1};
  G4int theAngularRep{-1};
  G4int nEnergy{-1};
  G4int theInterpolation{-1};

  G4InterpolationManager theManager;  // knows the interpolation between stores
  G4ParticleHPContAngularPar* theAngular;

  G4Cache<G4double> currentMeanEnergy;
  G4Cache<G4ParticleHPContAngularPar*> fCacheAngular;
  const G4ParticleDefinition* theProjectile;
};

#endif
