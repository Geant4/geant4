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

#ifndef G4ParticleHPInelasticBaseFS_h
#define G4ParticleHPInelasticBaseFS_h 1

#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4ParticleHPAngular.hh"
#include "G4ParticleHPDeExGammas.hh"
#include "G4ParticleHPEnAngCorrelation.hh"
#include "G4ParticleHPEnergyDistribution.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPPhotonDist.hh"
#include "G4ParticleHPNBodyPhaseSpace.hh"
#include "globals.hh"

class G4ParticleHPInelasticBaseFS : public G4ParticleHPFinalState
{
public:

  G4ParticleHPInelasticBaseFS();
  ~G4ParticleHPInelasticBaseFS() override;

  void Init(G4double A, G4double Z, G4int M, G4String& dirName, G4String& bit,
	    G4ParticleDefinition*) override;

  void BaseApply(const G4HadProjectile& theTrack, G4ParticleDefinition** theDefs, G4int nDef);

  void InitGammas(G4double AR, G4double ZR);

  G4HadFinalState* ApplyYourself(const G4HadProjectile& theTrack) override = 0;

  G4ParticleHPFinalState* New() override = 0;

  G4double GetXsec(G4double anEnergy) override
  {
    return std::max(0., theXsection->GetY(anEnergy));
  }

  G4ParticleHPVector* GetXsec() override { return theXsection; }

  G4ParticleHPInelasticBaseFS& operator=
  (const G4ParticleHPInelasticBaseFS& right) = delete;
  G4ParticleHPInelasticBaseFS(const G4ParticleHPInelasticBaseFS&) = delete;

protected:

  G4ParticleHPVector* theXsection;
  G4ParticleHPEnergyDistribution* theEnergyDistribution{nullptr};
  G4ParticleHPAngular* theAngularDistribution{nullptr};
  G4ParticleHPEnAngCorrelation* theEnergyAngData{nullptr};

  G4ParticleHPPhotonDist* theFinalStatePhotons{nullptr};
  G4double theNuclearMassDifference{0.0};
  G4double Qvalue{0.0};

  G4ParticleHPDeExGammas theGammas;
  G4String gammaPath{""};
};

#endif
