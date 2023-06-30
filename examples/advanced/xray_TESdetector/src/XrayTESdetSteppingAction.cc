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
/// \file XrayTESdetSteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XrayTESdetSteppingAction.hh"
#include "G4AnalysisManager.hh"

#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String clean_name (const G4String& vol_name)
{
  G4String cleaned_name = "";

  // Search if the vol_name contains one of the following part
  if (((std::string)vol_name).find("Bipxl") != std::string::npos)
    cleaned_name = "Bipxl";
  if (((std::string)vol_name).find("membrane") != std::string::npos)
    cleaned_name = "membrane";
  if (((std::string)vol_name).find("gridpiece") != std::string::npos)
    cleaned_name = "gridpiece";
  if (((std::string)vol_name).find("trapezoid") != std::string::npos)
    cleaned_name = "Mesh";
  if (cleaned_name == "") cleaned_name = vol_name;
  return cleaned_name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetSteppingAction::UserSteppingAction(const G4Step* step)
{
  // Define quantities to use in the
  G4RunManager *rm = G4RunManager::GetRunManager();
  auto analysisManager = G4AnalysisManager::Instance();

  G4int eventID = rm->GetCurrentEvent()->GetEventID();
  G4Track* track = step->GetTrack();
  G4StepPoint* pre_step_point = step->GetPreStepPoint();
  G4int parentID = track->GetParentID();
  G4int trackID = track->GetTrackID();
  const G4String& vol_name = track->GetVolume()->GetName();
  G4String mother_name = "";
  G4double init_kinetic_energy = 0;
  G4double kinetic_energy = track->GetKineticEnergy();
  G4double pre_kinetic_energy = pre_step_point->GetKineticEnergy();
  const G4String& particle_name = track->GetParticleDefinition()->GetParticleName();
  const G4VProcess* pre_step_proc = pre_step_point->GetProcessDefinedStep();
  const G4VProcess* post_step_proc = pre_step_point->GetProcessDefinedStep();

  bool proceed = false;
  G4String pre_step_name_s = "Undefined";
  G4String post_step_name_s = "Undefined";
  G4String creator_process_name = "Undefined";
  G4String init_creator_process = "Undefined";
  const G4String& pre_step_name = pre_step_name_s;
  const G4String& post_step_name = post_step_name_s;

  if (track->GetCreatorProcess() != nullptr)
  {
    creator_process_name = track->GetCreatorProcess()->GetProcessName();
  }

  if (pre_step_proc != nullptr)
  {
    pre_step_name_s = pre_step_proc->GetProcessName();
  }

  if (post_step_proc != nullptr)
  {
    post_step_name_s = post_step_proc->GetProcessName();
  }

  // Saves the first energy of the secondaries
  if (eventID != fPrev_eventID)
  {
    fInit_energy.clear();
    fCreator_proc.clear();
  }

  G4int current_step_number = -1;
  current_step_number = track->GetCurrentStepNumber();
  if (current_step_number != -1)
  {
    if (current_step_number == 1 && parentID == 0)
    {
      proceed = true;
      fInit_energy.insert(std::pair<G4int, G4double>(trackID, pre_kinetic_energy));
      fCreator_proc.insert(std::pair<G4int, G4String>(trackID, creator_process_name));
      pre_step_name_s = "InitStep";
    }
  }

  // Fill the map for the secondaries
  if (trackID > 1)
  {
    // First value for each trackIDs; this serves to keep it untile the eventID changes.
    if (fInit_energy.count((int)trackID) == 0)
    {
      fInit_energy.insert(std::pair<G4int, G4double>(trackID, pre_kinetic_energy));
      fCreator_proc.insert(std::pair<G4int, G4String>(trackID, creator_process_name));
    }
  }

  // Exclude the world volume, which has no mother volume
  if (((std::string)vol_name).find("expHall") == std::string::npos)
  {
    if (track->GetVolume() != nullptr)
    {
      if (track->GetVolume()->GetMotherLogical() != nullptr)
      {
        mother_name = track->GetVolume()->GetMotherLogical()->GetName();
      }
    }
    if (((std::string)mother_name).find("Detector_log") != std::string::npos || ((std::string)mother_name).find("element_log") != std::string::npos)
    {
      proceed = true;
    }
  }

  init_kinetic_energy = fInit_energy[(int)trackID];
  init_creator_process = fCreator_proc[trackID];

  G4String cleaned_name = "";
  cleaned_name = clean_name(vol_name);

  // Saves on tuple
  if (proceed)
  {
    G4double x = pre_step_point->GetPosition().x();
    G4double y = pre_step_point->GetPosition().y();
    G4double z = pre_step_point->GetPosition().z();
    G4ThreeVector direction = pre_step_point->GetMomentum();
    G4double theta = direction.theta();
    G4double phi = direction.phi();
    G4int pixel_number = track->GetVolume()->GetCopyNo();
    G4int def_replica_number = -999;                        // default value, to be discarded during the analysis
    G4int step_number = track->GetCurrentStepNumber();
    G4double step_energy_dep = step->GetTotalEnergyDeposit();

    // Fill tuple
    analysisManager->FillNtupleIColumn(0, eventID);
    analysisManager->FillNtupleSColumn(1, cleaned_name);
    analysisManager->FillNtupleIColumn(2, trackID);
    analysisManager->FillNtupleDColumn(3, x);
    analysisManager->FillNtupleDColumn(4, y);
    analysisManager->FillNtupleDColumn(5, z);
    analysisManager->FillNtupleDColumn(6, theta);
    analysisManager->FillNtupleDColumn(7, phi);
    analysisManager->FillNtupleIColumn(8, parentID);
    if (((std::string)vol_name).find("Bipxl") != std::string::npos)
    {
      analysisManager->FillNtupleIColumn(9, pixel_number);
    }
    else
    {
      analysisManager->FillNtupleIColumn(9, def_replica_number);
    }
    analysisManager->FillNtupleDColumn(10, step_energy_dep);
    analysisManager->FillNtupleIColumn(11, step_number);
    analysisManager->FillNtupleDColumn(12, init_kinetic_energy);
    analysisManager->FillNtupleDColumn(13, kinetic_energy);
    analysisManager->FillNtupleSColumn(14, particle_name);
    analysisManager->FillNtupleSColumn(15, pre_step_name);
    analysisManager->FillNtupleSColumn(16, post_step_name);
    analysisManager->FillNtupleSColumn(17, init_creator_process);
    analysisManager->AddNtupleRow(0);
  }
  fPrev_eventID = eventID;
  fPrev_trackID = trackID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
