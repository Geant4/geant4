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
/// \file radiobiology/src/Hit.cc
/// \brief Implementation of the RadioBio::Hit class

#include "Hit.hh"

#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

// for touchable handle
#include "G4TouchableHandle.hh"

#include <iomanip>

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<Hit>* RadioBioHitAllocator = 0;

Hit::Hit() : G4VHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Hit::Hit(const Hit& right) : G4VHit()
{
  fTrackID = right.fTrackID;
  fParticleType = right.fParticleType;
  fEkinMean = right.fEkinMean;
  fDeltaE = right.fDeltaE;
  fEinit = right.fEinit;
  fTrackLength = right.fTrackLength;
  fPhysicalVolume = right.fPhysicalVolume;
  fxIndex = right.fxIndex;
  fyIndex = right.fyIndex;
  fzIndex = right.fzIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Hit& Hit::operator=(const Hit& right)
{
  fTrackID = right.fTrackID;
  fParticleType = right.fParticleType;
  fEkinMean = right.fEkinMean;
  fDeltaE = right.fDeltaE;
  fEinit = right.fEinit;
  fTrackLength = right.fTrackLength;
  fPhysicalVolume = right.fPhysicalVolume;
  fxIndex = right.fxIndex;
  fyIndex = right.fyIndex;
  fzIndex = right.fzIndex;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int Hit::operator==(const Hit& right) const
{
  return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Hit::Draw()
{
  // Not implemented
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Hit::Print()
{
  // Not implemented
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Hit::SetVoxelIndexes(const G4TouchableHandle& TH)
{
  // Calculation of voxel number
  fxIndex = TH->GetReplicaNumber(2);
  fyIndex = TH->GetReplicaNumber(1);
  fzIndex = TH->GetReplicaNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio