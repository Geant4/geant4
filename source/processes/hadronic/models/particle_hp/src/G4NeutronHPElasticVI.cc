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
// Geant4 header : G4NeutronHPElasticVI
// Created:  15 October 2023
//
// Author  V.Ivanchenko
//

#include "G4NeutronHPElasticVI.hh"
#include "G4ParticleHPChannel.hh"
#include "G4ParticleHPElasticFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4Element.hh"
#include "G4Nucleus.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

G4bool G4NeutronHPElasticVI::fLock = false;
G4ParticleHPChannel* G4NeutronHPElasticVI::theElastic[] = {nullptr};

namespace
{
  G4Mutex theHPElastic = G4MUTEX_INITIALIZER;
}

G4NeutronHPElasticVI::G4NeutronHPElasticVI()
  : G4HadronicInteraction("NeutronHPElastic")
{
  SetMaxEnergy(20*CLHEP::MeV);
  fManagerHP = G4ParticleHPManager::GetInstance();
  if ( !fLock ) {
    fLock = true;
    fInitializer = true;
    for ( G4int i=0; i<ZMAXHPE; ++i ) {
      theElastic[i] = nullptr;
    }
  }
}

G4NeutronHPElasticVI::~G4NeutronHPElasticVI()
{
  if ( fInitializer ) {
    for ( G4int i=0; i<ZMAXHPE; ++i) {
      delete theElastic[i];
    }
  }
}

G4HadFinalState* G4NeutronHPElasticVI::ApplyYourself(const G4HadProjectile& aTrack,
                                                     G4Nucleus& aNucleus)
{
  G4HadFinalState* finalState = nullptr;
  G4int Z = aNucleus.GetZ_asInt();
  if ( Z >= ZMAXHPE || Z < 1 ) { return finalState; }

  G4int A = aNucleus.GetA_asInt();
  fManagerHP->OpenReactionWhiteBoard();
  fManagerHP->GetReactionWhiteBoard()->SetTargZ(Z);
  fManagerHP->GetReactionWhiteBoard()->SetTargA(A);

  G4ParticleHPChannel* mod = theElastic[Z];
  if ( nullptr == mod ) {
    InitialiseOnFly();
    if ( nullptr == mod ) { return finalState; }
  }

  // The boolean "true", as last argument, specifies to G4ParticleHPChannel::ApplyYourself
  // that it is an elastic channel: this is needed for the special DBRC treatment.
  finalState = mod->ApplyYourself(aTrack, -1, true);

  fManagerHP->CloseReactionWhiteBoard();
  return finalState;
}

const std::pair<G4double, G4double> G4NeutronHPElasticVI::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus
  return std::pair<G4double, G4double>(10.0 * perCent, 350.0 * CLHEP::GeV);
}

void G4NeutronHPElasticVI::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if ( fInitializer ) { Initialise(); }
}

void G4NeutronHPElasticVI::InitialiseOnFly()
{
  G4AutoLock l(&theHPElastic);
  Initialise();
  l.unlock();
}

void G4NeutronHPElasticVI::Initialise()
{
  G4ParticleHPElasticFS* theFS = nullptr;
  G4String dirName;
  for (auto const & elm : *(G4Element::GetElementTable())) {
    G4int Z = elm->GetZasInt();
    if ( 0 < Z && Z < ZMAXHPE && nullptr == theElastic[Z] ) {
      theElastic[Z] = new G4ParticleHPChannel();
      if ( nullptr == theFS ) {
	theFS = new G4ParticleHPElasticFS();
        dirName = fManagerHP->GetNeutronHPPath() + "/Elastic";
      }
      theElastic[Z]->Init(elm, dirName);
      theElastic[Z]->Register(theFS);
    }
  }
  delete theFS;
}

void G4NeutronHPElasticVI::ModelDescription(std::ostream& outFile) const
{
  outFile << "High Precision model based on Evaluated Nuclear Data Files"
          << " (ENDF) for elastic scattering of neutrons below 20MeV\n";
}
