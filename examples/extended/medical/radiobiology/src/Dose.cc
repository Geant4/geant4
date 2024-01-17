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
/// \file radiobiology/src/Dose.cc
/// \brief Implementation of the RadioBio::Dose class

#include "Dose.hh"

#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
#include "DoseAccumulable.hh"
#include "DoseMessenger.hh"
#include "VoxelizedSensitiveDetector.hh"

#include <fstream>

namespace RadioBio
{

#define width 15L

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Dose::Dose()
{
  // Default output filename
  fPath = "RadioBio_Dose.out";

  // Create the messenger
  fMessenger = new DoseMessenger(this);

  // Initialize the class
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Dose::~Dose()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Initialize()
{
  if (fVerboseLevel > 0) G4cout << "Dose::Initialize() called" << G4endl;

  G4int VoxelNumber = VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber();

  fEnDep.resize(VoxelNumber);
  fDose.resize(VoxelNumber);
  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Compute()
{
  // Skip Dose computation if calculation not enabled.
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "Dose::Compute() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  if (fCalculated) return;

  if (fVerboseLevel > 0) G4cout << "Dose::Compute() called" << G4endl;

  auto voxSensDet = VoxelizedSensitiveDetector::GetInstance();

  if (voxSensDet == nullptr)
    G4Exception("Dose::Compute", "VoxNotInit", FatalException,
                "Calculating dose without voxelized geometry pointer!");

  G4double voxelMass = voxSensDet->GetVoxelMass();

  for (G4int v = 0; v < voxSensDet->GetTotalVoxelNumber(); v++) {
    if (fVerboseLevel > 1) G4cout << "Calculating Dose for voxel number " << v << G4endl;

    // compute total Dose
    fDose[v] = fEnDep[v] / voxelMass;
  }

  fCalculated = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Store()
{
  // Skip Dose store if calculation not enabled.
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "Dose::Store() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  if (fSaved == true)
    G4Exception("Dose::Store", "DoseOverwrite", JustWarning,
                "Overwriting Dose file. For multiple runs, change filename.");

  Compute();

  std::ofstream ofs(fPath);

  if (ofs.is_open()) {
    ofs << "x_index" << std::setw(width) << "y_index" << std::setw(width) << "z_index"
        << std::setw(width) << "Dose (Gy)" << G4endl;

    auto voxSensDet = VoxelizedSensitiveDetector::GetInstance();

    for (G4int i = 0; i < voxSensDet->GetVoxelNumberAlongX(); i++)
      for (G4int j = 0; j < voxSensDet->GetVoxelNumberAlongY(); j++)
        for (G4int k = 0; k < voxSensDet->GetVoxelNumberAlongZ(); k++) {
          G4int v = voxSensDet->GetThisVoxelNumber(i, j, k);
          ofs << i << std::setw(width) << j << std::setw(width) << k << std::setw(width)
              << fDose[v] / gray << G4endl;
        }
  }
  if (fVerboseLevel > 0) {
    G4cout << "Dose: Dose written to " << fPath << G4endl;
  }
  fSaved = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::AddFromAccumulable(G4VAccumulable* GenAcc)
{
  DoseAccumulable* acc = (DoseAccumulable*)GenAcc;
  AddEnergyDeposit(acc->GetEnDeposit());
  fCalculated = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Reset()
{
  if (fVerboseLevel > 1) {
    G4cout << "Dose::Reset(): ";
  }
  fEnDep = 0.0;
  fDose = 0.0;
  fCalculated = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::PrintParameters()
{
  G4cout << "*******************************************" << G4endl
         << "****** Parameters of the class Dose *******" << G4endl
         << "*******************************************" << G4endl;
  PrintVirtualParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio