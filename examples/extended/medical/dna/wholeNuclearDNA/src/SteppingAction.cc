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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

// G4
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>
#include <G4Track.hh>

#include "Analysis.hh"
#include "CommandLineParser.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : G4UserSteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4double flagParticle = 0.;
  G4double flagProcess = 0.;
  G4double flagVolume = 0.;
  G4double x, y, z, xp, yp, zp;
  G4double dE;

  dE = step->GetTotalEnergyDeposit() / eV;

  const G4String& particleName = step->GetTrack()->GetDynamicParticle()
      ->GetDefinition()->GetParticleName();

  const G4String& processName =
      step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  const G4String& volumeName =
      step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if (particleName == "e-") flagParticle = 10;
  else if (particleName == "proton") flagParticle = 20;
  else if (particleName == "hydrogen") flagParticle = 30;
  else if (particleName == "alpha") flagParticle = 40;
  else if (particleName == "alpha+") flagParticle = 50;
  else if (particleName == "helium") flagParticle = 60;

  if (processName == "e-_G4DNAElastic") flagProcess = 11;
  else if (processName == "e-_G4DNAExcitation") flagProcess = 12;
  else if (processName == "e-_G4DNAIonisation") flagProcess = 13;
  else if (processName == "e-_G4DNAAttachment") flagProcess = 14;
  else if (processName == "e-_G4DNAVibExcitation") flagProcess = 15;
  else if (processName == "eCapture") flagProcess = 16;
//  if (step->GetPostStepPoint()->GetProcessDefinedStep()
// ->GetProcessName()=="msc")        flagProcess =17;

  else if (processName == "proton_G4DNAExcitation") flagProcess = 21;
  else if (processName == "proton_G4DNAIonisation") flagProcess = 22;
  else if (processName == "proton_G4DNAChargeDecrease") flagProcess = 23;

  else if (processName == "hydrogen_G4DNAExcitation") flagProcess = 31;
  else if (processName == "hydrogen_G4DNAIonisation") flagProcess = 32;
  else if (processName == "hydrogen_G4DNAChargeIncrease") flagProcess = 33;

  else if (processName == "alpha_G4DNAExcitation") flagProcess = 41;
  else if (processName == "alpha_G4DNAIonisation") flagProcess = 42;
  else if (processName == "alpha_G4DNAChargeDecrease") flagProcess = 43;

  else if (processName == "alpha+_G4DNAExcitation") flagProcess = 51;
  else if (processName == "alpha+_G4DNAIonisation") flagProcess = 52;
  else if (processName == "alpha+_G4DNAChargeDecrease") flagProcess = 53;
  else if (processName == "alpha+_G4DNAChargeIncrease") flagProcess = 54;

  else if (processName == "helium_G4DNAExcitation") flagProcess = 61;
  else if (processName == "helium_G4DNAIonisation") flagProcess = 62;
  else if (processName == "helium_G4DNAChargeIncrease") flagProcess = 63;

// if (step->GetPreStepPoint()->GetProcessDefinedStep()->
//  GetProcessName()=="hIoni")        flagProcess =24;
// if (step->GetPreStepPoint()->GetProcessDefinedStep()->
//  GetProcessName()=="eIoni")        flagProcess =18;

  if (volumeName == "physi sugar 2") flagVolume = 1;
  else if (volumeName == "physi sugar 4") flagVolume = 2;

  if (flagVolume != 0 && dE != 0)
  {

    x = step->GetPreStepPoint()->GetPosition().x() / nanometer;
    y = step->GetPreStepPoint()->GetPosition().y() / nanometer;
    z = step->GetPreStepPoint()->GetPosition().z() / nanometer;
    xp = step->GetPostStepPoint()->GetPosition().x() / nanometer;
    yp = step->GetPostStepPoint()->GetPosition().y() / nanometer;
    zp = step->GetPostStepPoint()->GetPosition().z() / nanometer;

    // The lines below could be put upper to gain time in the simulation
    // Added here for testing that all the retrieve information are 
    // correctly working
    CommandLineParser* parser = CommandLineParser::GetParser();
    Command* command(0);
    if((command = parser->GetCommandIfActive("-out"))==0) return;

    // get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleDColumn(0, flagParticle);
    analysisManager->FillNtupleDColumn(1, flagProcess);
    analysisManager->FillNtupleDColumn(2, flagVolume);
    analysisManager->FillNtupleDColumn(3, xp);
    analysisManager->FillNtupleDColumn(4, yp);
    analysisManager->FillNtupleDColumn(5, zp);
    analysisManager->FillNtupleDColumn(6, dE);
    analysisManager->FillNtupleDColumn(7,
                                       std::sqrt((x - xp) * (x - xp)
                                           + (y - yp) * (y - yp)
                                           + (z - zp) * (z - zp)));

    analysisManager->AddNtupleRow();
  }

}
