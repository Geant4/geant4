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
//---------------------------------------------------------------------------
//
// ClassName:   G4SpinDecayPhysics
//
// Author: 2015 P. Gumplinger
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4SpinDecayPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessTable.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

#include "G4EmBuilder.hh"
#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4SpinDecayPhysics);

G4SpinDecayPhysics::G4SpinDecayPhysics(G4int verb)
  : G4SpinDecayPhysics("SpinDecay", verb)
{
}

G4SpinDecayPhysics::G4SpinDecayPhysics(const G4String& name, G4int)
  : G4VPhysicsConstructor(name)
{
  G4PhysListUtil::InitialiseParameters();
}

G4SpinDecayPhysics::~G4SpinDecayPhysics()
{
}

void G4SpinDecayPhysics::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();

  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  
  G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
  MuonPlusDecayTable -> Insert(new
                         G4MuonDecayChannelWithSpin("mu+",0.986));
  MuonPlusDecayTable -> Insert(new
                         G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
  G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);

  G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
  MuonMinusDecayTable -> Insert(new
                          G4MuonDecayChannelWithSpin("mu-",0.986));
  MuonMinusDecayTable -> Insert(new
                          G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
  G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);
}

void G4SpinDecayPhysics::ConstructProcess()
{
  G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();

  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

  G4VProcess* decay;
  decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

  G4ProcessManager* fManager = G4MuonPlus::MuonPlus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(decayWithSpin);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
    fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
  }

  fManager = G4MuonMinus::MuonMinus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(decayWithSpin);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
    fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
  }

  G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

  decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

  fManager = G4PionPlus::PionPlus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }

  decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

  fManager = G4PionMinus::PionMinus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }
}
