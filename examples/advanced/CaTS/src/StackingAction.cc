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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file StackingAction.cc
/// \brief Implementation of the CaTS::StackingAction class

// Geant4 headers
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4GenericMessenger.hh"
// project headers
#include "StackingAction.hh"

StackingAction::StackingAction()
  : G4UserStackingAction()
{
  DefineCommands();
}

StackingAction::~StackingAction() { delete fMessenger; }
G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(
  const G4Track* aTrack)
{
  G4ClassificationOfNewTrack classification = fWaiting;
  if(aTrack->GetParentID() == 0)
    return classification;
  if(aTrack->GetDefinition()->GetParticleName() == "pi0")
  {
    if(fkillPi0)
    {
      classification = fKill;
    }
  }
  if(aTrack->GetDefinition()->GetParticleName() == "eta")
  {
    if(fkilleta)
    {
      classification = fKill;
    }
  }
  if(aTrack->GetDefinition()->GetParticleName() == "gamma")
  {
    if(aTrack->GetCreatorProcess()->GetProcessName() == "nCapture")
    {
      if(fkillGammafromnCapture)
      {
        classification = fKill;
      }
    }
  }
  return classification;
}

void StackingAction::Print()
{
  G4cout << "===================================================" << G4endl;
  G4cout << " StackingAction configuration:      " << G4endl;
  G4cout << " Kill Pi0s :                        " << fkillPi0 << G4endl;
  G4cout << " Kill etas :                        " << fkilleta << G4endl;
  G4cout << " Kill Gammas from neutron Capture:  " << fkillGammafromnCapture
         << G4endl;
  G4cout << "===================================================" << G4endl;
}

void StackingAction::DefineCommands()
{
  fMessenger       = new G4GenericMessenger(this, "/CaTS/StackingAction/",
                                      "select particles to kill");
  auto& killPi0Cmd = fMessenger->DeclareProperty("killPi0", fkillPi0);
  killPi0Cmd.SetGuidance("kill Pi0 (true/false)");
  killPi0Cmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  killPi0Cmd.SetDefaultValue("false");
  auto& killetaCmd = fMessenger->DeclareProperty("killeta", fkilleta);
  killetaCmd.SetGuidance("kill eta (true/false)");
  killetaCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  killetaCmd.SetDefaultValue("false");
  auto& killGammafromnCaptureCmd = fMessenger->DeclareProperty(
    "killGammafromnCapture", fkillGammafromnCapture);
  killGammafromnCaptureCmd.SetGuidance("kill GammafromnCapture (true/false)");
  killGammafromnCaptureCmd.SetStates(G4State_PreInit, G4State_Init,
                                     G4State_Idle);
  killGammafromnCaptureCmd.SetDefaultValue("false");
  fMessenger->DeclareMethod("Print", &StackingAction::Print)
    .SetGuidance("Print StackingAction configuration")
    .SetStates(G4State_PreInit, G4State_Idle);
}
