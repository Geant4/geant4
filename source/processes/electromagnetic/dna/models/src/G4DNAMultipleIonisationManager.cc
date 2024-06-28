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
// G4DNAMultipleIonisationManager.cc
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//

#include "G4DNAMultipleIonisationManager.hh"
#include "G4FindDataDir.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "G4VITTrackHolder.hh"
#include "G4H2O.hh"
#include "G4DNAChemistryManager.hh"

//------------------------------------------------------------------------------
void G4DNAMultipleIonisationManager::CreateMultipleIonisedWaterMolecule(
  MultipleIonisedModification mod, G4int* shell_level,
  const G4Track* incoming_track)
{
  if (!G4DNAChemistryManager::IsActivated()) { return; }

  G4int num_shells{0};
  switch (mod) {
    case eDoubleIonisedMolecule:
      num_shells = 2;
      break;
    case eTripleIonisedMolecule:
      num_shells = 3;
      break;
    case eQuadrupleIonisedMolecule:
      num_shells = 4;
      break;
    default: // never happen
      return;
  }

  auto* H2O = new G4Molecule(G4H2O::Definition());
  for (G4int i = 0; i < num_shells; i++) {
    H2O->IonizeMolecule(4 - shell_level[i]);
  }

  constexpr G4double kT0 = 1.0 * picosecond;
  auto* H2O_track = H2O->BuildTrack(kT0, incoming_track->GetPosition());
  H2O_track->SetParentID(incoming_track->GetTrackID());
  H2O_track->SetTrackStatus(fStopButAlive);
  H2O_track->SetKineticEnergy(0.0);

  G4VITTrackHolder::Instance()->Push(H2O_track);
}

//------------------------------------------------------------------------------
G4bool G4DNAMultipleIonisationManager::CheckShellEnergy(
    MultipleIonisedModification mod, G4double* shell_energy)
{
  G4int num_shells{0};
  switch (mod) {
    case eDoubleIonisedMolecule:
      num_shells = 2;
      break;
    case eTripleIonisedMolecule:
      num_shells = 3;
      break;
    case eQuadrupleIonisedMolecule:
      num_shells = 4;
      break;
    default: // never happen
      break;
  }

  G4bool stop_process{false};
  for (int i = 0; i < num_shells; i++) {
    if (shell_energy[i] < 0.0) {
      stop_process = true;
      break;
    }
  }

  return stop_process;
}

//------------------------------------------------------------------------------
void G4DNAMultipleIonisationManager::LoadAlphaParam(
  const G4String& filepath, G4double Z, G4double A)
{
  const char* path = G4FindDataDir("G4LEDATA");
  if (path == nullptr) {
    G4Exception("G4DNAMultipleIonisationManager::LoadAlphaParam","em0006",
                FatalException,"G4LEDATA environment variable not set.");
  }

  std::stringstream fullpath;
  fullpath << path << "/" << filepath;

  std::fstream fin(fullpath.str());

  G4double e, a;
  std::string line = "";
  while (getline(fin, line)) {
    std::stringstream ss;
    ss << line;
    ss >> e >> a;
    Eion_.push_back(e * Z * A * MeV);
    alpha_.push_back(a);
  }

  num_node_ = (G4int)Eion_.size();
  fin.close();
}

//------------------------------------------------------------------------------
G4double G4DNAMultipleIonisationManager::GetAlphaParam(G4double energy)
{
  auto find_lower_bound = [this](G4double e) {
    auto low = 0;
    auto upp = num_node_ - 1;
    if (e < Eion_[0]) { return low; }
    while (low <= upp) {
      const auto mid = static_cast<int>((low + upp) * 0.5);
      if (e < Eion_[mid]) { upp = mid - 1; }
      else { low = mid + 1; }
      if (upp < 0) { upp = 0; }
    }
    return upp;
  };

  auto interp_log_log = [this](G4int bin1, G4double e) {
    if (e < Eion_[0]) { return alpha_[0]; }
    const auto num_bin = num_node_ - 1;
    const auto bin2 = bin1 + 1;
    G4double value{0.0};
    if (bin2 <= num_bin) {
      auto log10_e  = std::log10(e);
      auto log10_e1 = Eion_[bin1];
      auto log10_e2 = Eion_[bin2];
      auto log10_a1 = alpha_[bin1];
      auto log10_a2 = alpha_[bin2];
      if (log10_a1 != 0.0 && log10_a2 != 0.0) {
        log10_e1 = std::log10(log10_e1);
        log10_e2 = std::log10(log10_e2);
        log10_a1 = std::log10(log10_a1);
        log10_a2 = std::log10(log10_a2);
        value = log10_a1 + (log10_a2 - log10_a1)
                  * (log10_e - log10_e1) / (log10_e2 - log10_e1);
        value = std::pow(10.0, value);
      }
    } else {
      value = alpha_[num_bin];
    }
    return value;
  };

  const auto bin1 = find_lower_bound(energy);
  return interp_log_log(bin1, energy);
}
