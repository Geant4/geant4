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
// Geant4 class : G4NeutronHPInelasticVI
// Created:  15 October 2023
//
// Author  V.Ivanchenko
//

#include "G4NeutronHPInelasticVI.hh"

#include "G4HadronicParameters.hh"
#include "G4ParticleHPManager.hh"

#include "G4ParticleHP2AInelasticFS.hh"
#include "G4ParticleHP2N2AInelasticFS.hh"
#include "G4ParticleHP2NAInelasticFS.hh"
#include "G4ParticleHP2NDInelasticFS.hh"
#include "G4ParticleHP2NInelasticFS.hh"
#include "G4ParticleHP2NPInelasticFS.hh"
#include "G4ParticleHP2PInelasticFS.hh"
#include "G4ParticleHP3AInelasticFS.hh"
#include "G4ParticleHP3NAInelasticFS.hh"
#include "G4ParticleHP3NInelasticFS.hh"
#include "G4ParticleHP3NPInelasticFS.hh"
#include "G4ParticleHP4NInelasticFS.hh"
#include "G4ParticleHPAInelasticFS.hh"
#include "G4ParticleHPD2AInelasticFS.hh"
#include "G4ParticleHPDAInelasticFS.hh"
#include "G4ParticleHPDInelasticFS.hh"
#include "G4ParticleHPHe3InelasticFS.hh"
#include "G4ParticleHPN2AInelasticFS.hh"
#include "G4ParticleHPN2PInelasticFS.hh"
#include "G4ParticleHPN3AInelasticFS.hh"
#include "G4ParticleHPNAInelasticFS.hh"
#include "G4ParticleHPND2AInelasticFS.hh"
#include "G4ParticleHPNDInelasticFS.hh"
#include "G4ParticleHPNHe3InelasticFS.hh"
#include "G4ParticleHPNInelasticFS.hh"
#include "G4ParticleHPNPAInelasticFS.hh"
#include "G4ParticleHPNPInelasticFS.hh"
#include "G4ParticleHPNT2AInelasticFS.hh"
#include "G4ParticleHPNTInelasticFS.hh"
#include "G4ParticleHPNXInelasticFS.hh"
#include "G4ParticleHPPAInelasticFS.hh"
#include "G4ParticleHPPDInelasticFS.hh"
#include "G4ParticleHPPInelasticFS.hh"
#include "G4ParticleHPPTInelasticFS.hh"
#include "G4ParticleHPT2AInelasticFS.hh"
#include "G4ParticleHPTInelasticFS.hh"
#include "G4ParticleHPThermalBoost.hh"

#include "G4Neutron.hh"
#include "G4Nucleus.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

G4bool G4NeutronHPInelasticVI::fLock = false;
G4ParticleHPChannelList* G4NeutronHPInelasticVI::theChannels[] = {nullptr};

namespace
{
  G4Mutex theHPInelastic = G4MUTEX_INITIALIZER;
}

G4NeutronHPInelasticVI::G4NeutronHPInelasticVI()
  : G4HadronicInteraction("NeutronHPInelastic")
{
  SetMaxEnergy(20*CLHEP::MeV);
  fManagerHP = G4ParticleHPManager::GetInstance();
  if ( !fLock ) {
    fLock = true;
    fInitializer = true;
    for ( G4int i=0; i<ZMAXHPI; ++i ) {
      theChannels[i] = nullptr;
    }
  }
}

G4NeutronHPInelasticVI::~G4NeutronHPInelasticVI()
{
  if ( fInitializer ) {
    for ( G4int i=0; i<ZMAXHPI; ++i ) {
      delete theChannels[i];
    }
  }
}

G4HadFinalState* G4NeutronHPInelasticVI::ApplyYourself(const G4HadProjectile& aTrack,
                                                       G4Nucleus& aNucleus)
{
  G4HadFinalState* finalState = nullptr;
  G4int Z = aNucleus.GetZ_asInt();
  if ( Z >= ZMAXHPI || Z < 1 ) { return finalState; }

  G4int A = aNucleus.GetA_asInt();
  fManagerHP->OpenReactionWhiteBoard();
  fManagerHP->GetReactionWhiteBoard()->SetTargZ(Z);
  fManagerHP->GetReactionWhiteBoard()->SetTargA(A);

  G4ParticleHPChannelList* clist = theChannels[Z];
  if ( nullptr == clist ) {
    InitialiseOnFly();
    if ( nullptr == clist ) { return finalState; }
  }

  for (auto const & elm : *(G4Element::GetElementTable())) {
    if ( Z == elm->GetZasInt() ) {
      finalState = clist->ApplyYourself(elm, aTrack);
      break;
    }
  }

  G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
  return finalState;
}

const std::pair<G4double, G4double> G4NeutronHPInelasticVI::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus
  return std::pair<G4double, G4double>(10.0 * perCent, 350.0 * CLHEP::GeV);
}

void G4NeutronHPInelasticVI::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if ( fInitializer ) {
    Initialise();
    fManagerHP->DumpSetting();
  }
}

void G4NeutronHPInelasticVI::InitialiseOnFly()
{
  G4AutoLock l(&theHPInelastic);
  Initialise();
  l.unlock();
}

void G4NeutronHPInelasticVI::Initialise()
{
  G4String dirName;
  G4ParticleDefinition* part = nullptr;

  for (auto const & elm : *(G4Element::GetElementTable())) {
    G4int Z = elm->GetZasInt();
    if ( 0 < Z && Z < ZMAXHPI && nullptr == theChannels[Z] ) {
      if ( nullptr == part ) {
        part = G4Neutron::Neutron();
        dirName = fManagerHP->GetNeutronHPPath() + "/Inelastic";
      }
      auto clist = new G4ParticleHPChannelList(36, part);
      theChannels[Z] = clist;
      clist->Init(elm, dirName, part);
      clist->Register(new G4ParticleHPNInelasticFS, "F01/");  // has
      clist->Register(new G4ParticleHPNXInelasticFS, "F02/");
      clist->Register(new G4ParticleHP2NDInelasticFS, "F03/");
      clist->Register(new G4ParticleHP2NInelasticFS, "F04/");  // has, E Done
      clist->Register(new G4ParticleHP3NInelasticFS, "F05/");  // has, E Done
      clist->Register(new G4ParticleHPNAInelasticFS, "F06/");
      clist->Register(new G4ParticleHPN3AInelasticFS, "F07/");
      clist->Register(new G4ParticleHP2NAInelasticFS, "F08/");
      clist->Register(new G4ParticleHP3NAInelasticFS, "F09/");
      clist->Register(new G4ParticleHPNPInelasticFS, "F10/");
      clist->Register(new G4ParticleHPN2AInelasticFS, "F11/");
      clist->Register(new G4ParticleHP2N2AInelasticFS, "F12/");
      clist->Register(new G4ParticleHPNDInelasticFS, "F13/");
      clist->Register(new G4ParticleHPNTInelasticFS, "F14/");
      clist->Register(new G4ParticleHPNHe3InelasticFS, "F15/");
      clist->Register(new G4ParticleHPND2AInelasticFS, "F16/");
      clist->Register(new G4ParticleHPNT2AInelasticFS, "F17/");
      clist->Register(new G4ParticleHP4NInelasticFS, "F18/");  // has, E Done
      clist->Register(new G4ParticleHP2NPInelasticFS, "F19/");
      clist->Register(new G4ParticleHP3NPInelasticFS, "F20/");
      clist->Register(new G4ParticleHPN2PInelasticFS, "F21/");
      clist->Register(new G4ParticleHPNPAInelasticFS, "F22/");
      clist->Register(new G4ParticleHPPInelasticFS, "F23/");
      clist->Register(new G4ParticleHPDInelasticFS, "F24/");
      clist->Register(new G4ParticleHPTInelasticFS, "F25/");
      clist->Register(new G4ParticleHPHe3InelasticFS, "F26/");
      clist->Register(new G4ParticleHPAInelasticFS, "F27/");
      clist->Register(new G4ParticleHP2AInelasticFS, "F28/");
      clist->Register(new G4ParticleHP3AInelasticFS, "F29/");
      clist->Register(new G4ParticleHP2PInelasticFS, "F30/");
      clist->Register(new G4ParticleHPPAInelasticFS, "F31/");
      clist->Register(new G4ParticleHPD2AInelasticFS, "F32/");
      clist->Register(new G4ParticleHPT2AInelasticFS, "F33/");
      clist->Register(new G4ParticleHPPDInelasticFS, "F34/");
      clist->Register(new G4ParticleHPPTInelasticFS, "F35/");
      clist->Register(new G4ParticleHPDAInelasticFS, "F36/");
#ifdef G4VERBOSE
      if (fManagerHP->GetVerboseLevel() > 1) {
	G4cout << "G4NeutronHP::InelasticVI for " 
	       << part->GetParticleName() << " off " 
	       << elm->GetName() << G4endl;
      }
#endif
    }
  }
}

void G4NeutronHPInelasticVI::ModelDescription(std::ostream& outFile) const
{
  outFile << "High Precision (HP) model for inelastic reaction of "
	  << " neutrons below 20MeV\n";
}
