// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LeptonConstructor.cc,v 1.2 1999-12-15 14:51:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementatLepton file 
//

#include "G4LeptonConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
// Leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
G4LeptonConstructor::G4LeptonConstructor()
{
}

G4LeptonConstructor::~G4LeptonConstructor()
{
}


void G4LeptonConstructor::ConstructParticle()
{
  ConstructELeptons();
  ConstructMuLeptons();
  ConstructTauLeptons();
}

void G4LeptonConstructor::ConstructELeptons()
{
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

void G4LeptonConstructor::ConstructMuLeptons()
{
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void G4LeptonConstructor::ConstructTauLeptons()
{
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
}
