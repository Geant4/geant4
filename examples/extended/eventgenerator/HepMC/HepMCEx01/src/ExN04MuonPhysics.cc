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
// $Id: ExN04MuonPhysics.cc,v 1.4 2003/12/03 11:03:06 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

#include "ExN04MuonPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   


ExN04MuonPhysics::ExN04MuonPhysics(const G4String& name)
                   :  G4VPhysicsConstructor(name)
{
}

ExN04MuonPhysics::~ExN04MuonPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

void ExN04MuonPhysics::ConstructParticle()
{
  // Mu
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // Tau
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();

}


#include "G4ProcessManager.hh"

void ExN04MuonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  // Muon Plus Physics
  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
   
  pManager->AddProcess(&fMuPlusMultipleScattering,-1,  1, 1);
  pManager->AddProcess(&fMuPlusIonisation,        -1,  2, 2);
  pManager->AddProcess(&fMuPlusBremsstrahlung,    -1,  3, 3);
  pManager->AddProcess(&fMuPlusPairProduction,    -1,  4, 4);

  // Muon Minus Physics
  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
   
  pManager->AddProcess(&fMuMinusMultipleScattering,-1,  1, 1);
  pManager->AddProcess(&fMuMinusIonisation,        -1,  2, 2);
  pManager->AddProcess(&fMuMinusBremsstrahlung,    -1,  3, 3);
  pManager->AddProcess(&fMuMinusPairProduction,    -1,  4, 4);


  pManager->AddRestProcess(&fMuMinusCaptureAtRest);

  // Tau Plus Physics
  pManager = G4TauPlus::TauPlus()->GetProcessManager();
   
  pManager->AddProcess(&fTauPlusMultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fTauPlusIonisation,         -1, 2, 2);

  // Tau Minus Physics
  pManager = G4TauMinus::TauMinus()->GetProcessManager();

  pManager->AddProcess(&fTauMinusMultipleScattering, -1, 1, 1);
  pManager->AddProcess(&fTauMinusIonisation,         -1, 2, 2);

}



