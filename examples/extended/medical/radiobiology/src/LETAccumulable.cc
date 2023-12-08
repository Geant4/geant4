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
/// \file radiobiology/src/LETAccumulable.cc
/// \brief Implementation of the RadioBio::LETAccumulable class

#include "LETAccumulable.hh"

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

LETAccumulable::LETAccumulable() : VRadiobiologicalAccumulable("LET") {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LETAccumulable::Merge(const G4VAccumulable& rhs)
{
  if (GetVerboseLevel() > 0) {
    G4cout << "LETAccumulable::Merge()" << G4endl;
  }
  const LETAccumulable& other = dynamic_cast<const LETAccumulable&>(rhs);

  // Merges standard counters
  fTotalLETT += other.GetTotalLETT();
  fDTotalLETT += other.GetDTotalLETT();
  fDTotalLETD += other.GetDTotalLETD();
  fTotalLETD += other.GetTotalLETD();

  // Merges ion counters
  auto rhsIonLetStore = other.GetIonLetStore();

  // Loop over rhs ions
  for (unsigned int l = 0; l < rhsIonLetStore.size(); l++) {
    G4int PDGencoding = rhsIonLetStore[l].GetPDGencoding();
    size_t q;
    // Loop over lhs ions to find the one
    for (q = 0; q < fIonLetStore.size(); q++) {
      // Check if primary status is the same
      if (fIonLetStore[q].GetPDGencoding() == PDGencoding) {
        if (rhsIonLetStore[l].IsPrimary() == fIonLetStore[q].IsPrimary()) break;
      }
    }

    if (q == fIonLetStore.size())  // Just copy the IonLet element in lhs store
      fIonLetStore.push_back(rhsIonLetStore[l]);
    else  // Merges rhs data with lhs ones
      fIonLetStore[q].Merge(&(rhsIonLetStore[l]));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LETAccumulable::Reset()
{
  if (GetVerboseLevel() > 0) {
    G4cout << "LETAccumulable::Reset()" << G4endl;
  }
  if (fInitialized) {
    fTotalLETT = 0.0;
    fDTotalLETT = 0.0;
    fDTotalLETD = 0.0;
    fTotalLETD = 0.0;
    fIonLetStore.clear();
  }
  else {
    Initialize();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// To accumulate given the hit
void LETAccumulable::Accumulate(Hit* hit)
{
  if (GetVerboseLevel() > 1) {
    G4cout << "LETAccumulable::Accumulate()" << G4endl;
  }

  // No need to accumulate if no energy deposit or track length is zero.
  if (hit->GetDeltaE() == 0. || hit->GetTrackLength() == 0.) return;

  // Atomic number
  G4int Z = hit->GetPartType()->GetAtomicNumber();
  if (Z < 1) return;  // Calculate only protons and ions

  G4int PDGencoding = hit->GetPartType()->GetPDGEncoding();

  if (PDGencoding == 22 || PDGencoding == 11) return;  // do not accumulate for electrons or gammas

  // Get simple Particle data Group ID without excitation level
  PDGencoding -= PDGencoding % 10;

  // Get voxel number
  G4int xIndex = hit->GetXindex();
  G4int yIndex = hit->GetYindex();
  ;
  G4int zIndex = hit->GetZindex();
  ;

  G4int voxel =
    VoxelizedSensitiveDetector::GetInstance()->GetThisVoxelNumber(xIndex, yIndex, zIndex);

  // Get mean kinetic energy
  G4double ekinMean = hit->GetEkinMean();

  // Get particle definition
  const G4ParticleDefinition* particleDef = hit->GetPartType();

  // Get material
  G4Material* mat = hit->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();

  // ICRU stopping power calculation
  G4EmCalculator emCal;
  // Use the mean kinetic energy of ions in a step to calculate ICRU stopping power
  G4double Lsn = emCal.ComputeElectronicDEDX(ekinMean, particleDef, mat);
  // G4cout << "ERASE ME. Lsn = " << Lsn << G4endl;
  // G4cout << "ERASE ME. ekinMean = " << ekinMean << G4endl;

  // Get deposited energy
  G4double DE = hit->GetDeltaE();

  // Get deposited energy considering secondary electrons
  G4double DEElectrons = hit->GetElectronEnergy();

  // Get track length
  G4double DX = hit->GetTrackLength();

  // Get track ID
  G4int trackID = hit->GetTrackID();

  // Total LET calculation...
  // Total dose LET Numerator, including secondary electrons energy deposit
  fTotalLETD[voxel] += (DE + DEElectrons) * Lsn;
  // Total dose LET Denominator, including secondary electrons energy deposit
  fDTotalLETD[voxel] += DE + DEElectrons;
  // Total track LET Numerator
  fTotalLETT[voxel] += DX * Lsn;
  // Total track LET Denominator
  fDTotalLETT[voxel] += DX;

  size_t l;
  for (l = 0; l < fIonLetStore.size(); l++) {
    // Judge species of the current particle and store it
    if (fIonLetStore[l].GetPDGencoding() == PDGencoding)
      if (((trackID == 1) && (fIonLetStore[l].IsPrimary()))
          || ((trackID != 1) && (!fIonLetStore[l].IsPrimary())))
        break;
  }

  // Just another type of ion/particle for our store...
  if (l == fIonLetStore.size()) {
    // Mass number
    G4int A = particleDef->GetAtomicMass();

    // Particle name
    G4String fullName = particleDef->GetParticleName();
    G4String name = fullName.substr(0, fullName.find("["));  // Cut excitation energy [x.y]

    IonLet ion(trackID, PDGencoding, fullName, name, Z, A,
               VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber());
    fIonLetStore.push_back(ion);
  }

  // Calculate ions LET and store them
  fIonLetStore[l].Update(voxel, DE, DEElectrons, Lsn, DX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int LETAccumulable::GetVerboseLevel() const
{
  // return same level of LET class
  return Manager::GetInstance()->GetQuantity("LET")->GetVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LETAccumulable::Initialize()
{
  if (GetVerboseLevel() > 0) {
    G4cout << "LETAccumulable::Initialize(): " << G4endl;
  }

  G4int voxNumber = VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber();

  fTotalLETT = array_type(0.0, voxNumber);
  fTotalLETD = array_type(0.0, voxNumber);
  fDTotalLETT = array_type(0.0, voxNumber);
  fDTotalLETD = array_type(0.0, voxNumber);
  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio