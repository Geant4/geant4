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
// $Id: ExN04MuonPhysics.cc,v 1.5 2006/06/29 17:06:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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



