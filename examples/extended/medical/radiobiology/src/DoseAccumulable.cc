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
/// \file radiobiology/src/DoseAccumulable.cc
/// \brief Implementation of the RadioBio::DoseAccumulable class

#include "DoseAccumulable.hh"

#include "G4EmCalculator.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"

#include "Hit.hh"
#include "Manager.hh"
#include "VoxelizedSensitiveDetector.hh"
#include <G4SystemOfUnits.hh>

#include <tuple>

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DoseAccumulable::DoseAccumulable() : VRadiobiologicalAccumulable("Dose") {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Merge(const G4VAccumulable& rhs)
{
  if (GetVerboseLevel() > 1) {
    G4cout << "DoseAccumulable::Merge()" << G4endl;
  }
  const DoseAccumulable& other = dynamic_cast<const DoseAccumulable&>(rhs);

  // Merges the counter
  fEnDep += other.GetEnDeposit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Reset()
{
  if (GetVerboseLevel() > 0) {
    G4cout << "DoseAccumulable::Reset()" << G4endl;
  }
  if (fInitialized) {
    fEnDep = 0.0;
  }
  else {
    Initialize();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// To accumulate given the hit
void DoseAccumulable::Accumulate(Hit* hit)
{
  // Calculation done only if energy is deposited
  if (hit->GetDeltaE() <= 0.) return;

  if (GetVerboseLevel() > 1) {
    G4cout << "DoseAccumulable::Accumulate()" << G4endl;
  }

  // Get voxel number
  G4int xIndex = hit->GetXindex();
  G4int yIndex = hit->GetYindex();
  G4int zIndex = hit->GetZindex();

  G4int voxel =
    VoxelizedSensitiveDetector::GetInstance()->GetThisVoxelNumber(xIndex, yIndex, zIndex);

  // Get deposited energy
  G4double DE = hit->GetDeltaE();

  // Update total dose
  fEnDep[voxel] += DE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DoseAccumulable::GetVerboseLevel() const
{
  // Returns same verbosity of Dose class
  return Manager::GetInstance()->GetQuantity("Dose")->GetVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Initialize()
{
  if (GetVerboseLevel() > 0) {
    G4cout << "DoseAccumulable::Initialize(): " << G4endl;
  }

  fEnDep = array_type(0.0, VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber());

  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio