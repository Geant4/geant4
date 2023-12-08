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
/// \file radiobiology/src/VoxelizedSensitiveDetector.cc
/// \brief Implementation of the RadioBio::VoxelizedSensitiveDetector class

#include "VoxelizedSensitiveDetector.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

#include "DetectorConstruction.hh"
#include "SD.hh"
#include "VoxelizedSensitiveDetectorMessenger.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector* VoxelizedSensitiveDetector::fInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector* VoxelizedSensitiveDetector::CreateInstance(DetectorConstruction* det,
                                                                       double xWidth, double yWidth,
                                                                       double zWidth)
{
  if (fInstance) {
    delete fInstance;
    G4Exception("VoxelizedSensitiveDetector::createInstance", "RecreatingVoxelization",
                FatalException, "Creating another, new, instance of VoxelizedSensitiveDetector");
  }
  fInstance = new VoxelizedSensitiveDetector(det, xWidth, yWidth, zWidth);
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector* VoxelizedSensitiveDetector::GetInstance()
{
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector::VoxelizedSensitiveDetector(DetectorConstruction* det, double xWidth,
                                                       double yWidth, double zWidth)
  : fDetector(det), fVoxelWidthX(xWidth), fVoxelWidthY(yWidth), fVoxelWidthZ(zWidth)
{
  fVoxelizedSensitiveDetectorMessenger = new VoxelizedSensitiveDetectorMessenger(this);
  UpdateVoxelVolume();
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector::~VoxelizedSensitiveDetector()
{
  delete fVoxelizedSensitiveDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::UpdateVoxelVolume()
{
  fVoxelVolume = fVoxelWidthX * fVoxelWidthY * fVoxelWidthZ;
  fVoxelDensity = fDetector->GetMaterial()->GetDensity();
  fVoxelMass = fVoxelVolume * fVoxelDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidth(G4ThreeVector voxWidth)
{
  fVoxelWidthX = voxWidth.getX();
  fVoxelWidthY = voxWidth.getY();
  fVoxelWidthZ = voxWidth.getZ();
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidthX(G4double voxWidthX)
{
  if (fVoxelWidthX == voxWidthX) return;
  fVoxelWidthX = voxWidthX;
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidthY(G4double voxWidthY)
{
  if (fVoxelWidthY == voxWidthY) return;
  fVoxelWidthY = voxWidthY;
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidthZ(G4double voxWidthZ)
{
  if (fVoxelWidthZ == voxWidthZ) return;
  fVoxelWidthZ = voxWidthZ;
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Calculte number of voxel approximating for an integer number of voxels
// Then recalculates voxels size according to approximations
void VoxelizedSensitiveDetector::CalculateVoxelNumber()
{
  fVoxelNumberAlongX = G4int(fDetector->GetSizeX() / fVoxelWidthX);
  fVoxelWidthX = fDetector->GetSizeX() / G4double(fVoxelNumberAlongX);

  fVoxelNumberAlongY = G4int(fDetector->GetSizeY() / fVoxelWidthY);
  fVoxelWidthY = fDetector->GetSizeY() / G4double(fVoxelNumberAlongY);

  fVoxelNumberAlongZ = G4int(fDetector->GetSizeZ() / fVoxelWidthZ);
  fVoxelWidthZ = fDetector->GetSizeZ() / G4double(fVoxelNumberAlongZ);

  if (fVoxelNumberAlongY % 2 == 0)
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "VoxelNumberYEven", JustWarning,
                "Trying to voxelize with an even number of voxels along the Y axis."
                "Please select an odd number to prevent from warnings due to tracking");

  if (fVoxelNumberAlongZ % 2 == 0)
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "VoxelNumberZEven", JustWarning,
                "Trying to voxelize with an even number of voxels along the Z axis."
                "Please select an odd number to prevent from warnings due to tracking");

  fTotalVoxelNumber = fVoxelNumberAlongX * fVoxelNumberAlongY * fVoxelNumberAlongZ;

  UpdateVoxelVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructXDivision()
{
  if (fWorldLogical == nullptr)
    G4Exception("VoxelizedSensitiveDetector::ConstructXDivision", "WorldNotInit", FatalException,
                "Voxelizing without having a pointer to world logical volume!");

  if (!fDetector)
    G4Exception("VoxelizedSensitiveDetector::ConstructXDivision", "DetConstInit", FatalException,
                "Voxelizing without having a pointer to DetectorConstruction!");

  fVoxelizedDetectorXDivision = new G4Box("VoxelizedDetectorXDivision", fVoxelWidthX / 2,
                                          fDetector->GetSizeY() / 2, fDetector->GetSizeZ() / 2);

  fVoxelizedDetectorXDivisionLog =
    new G4LogicalVolume(fVoxelizedDetectorXDivision, fWorldLogical->GetMaterial(),
                        "VoxelizedDetectorXDivisionLog", 0, 0, 0);

  fVoxelizedDetectorXDivisionPhys =
    new G4PVReplica("VoxelizedDetectorXDivisionPhys", fVoxelizedDetectorXDivisionLog, fWorldLogical,
                    kXAxis, fVoxelNumberAlongX, fVoxelWidthX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructYDivision()
{
  fVoxelizedDetectorYDivision = new G4Box("VoxelizedDetectorYDivision", fVoxelWidthX / 2,
                                          fVoxelWidthY / 2, fDetector->GetSizeZ() / 2);

  fVoxelizedDetectorYDivisionLog =
    new G4LogicalVolume(fVoxelizedDetectorYDivision, fWorldLogical->GetMaterial(),
                        "VoxelizedDetectorYDivisionLog", 0, 0, 0);

  fVoxelizedDetectorYDivisionPhys =
    new G4PVReplica("VoxelizedDetectorYDivisionPhys", fVoxelizedDetectorYDivisionLog,
                    fVoxelizedDetectorXDivisionLog, kYAxis, fVoxelNumberAlongY, fVoxelWidthY);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructZDivision()
{
  fVoxelizedDetectorZDivision =
    new G4Box("VoxelizedDetectorZDivision", fVoxelWidthX / 2, fVoxelWidthY / 2, fVoxelWidthZ / 2);

  fVoxelizedDetectorZDivisionLog =
    new G4LogicalVolume(fVoxelizedDetectorZDivision, fWorldLogical->GetMaterial(),
                        "VoxelizedDetectorZDivisionLog", 0, 0, 0);

  fVoxelizedDetectorZDivisionPhys =
    new G4PVReplica("VoxelizedDetectorZDivisionPhys", fVoxelizedDetectorZDivisionLog,
                    fVoxelizedDetectorYDivisionPhys, kZAxis, fVoxelNumberAlongZ, fVoxelWidthZ);

  fSensitiveLogicalVolume = fVoxelizedDetectorZDivisionLog;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// First voxelize along X, then Y, then Z
G4bool VoxelizedSensitiveDetector::ConstructVoxelizedDetector()
{
  // Creating X division
  ConstructXDivision();

  // Creating Y division
  ConstructYDivision();

  // Creating Z division
  ConstructZDivision();

  // Set last, smallest volumes as sensitive
  fSensitiveLogicalVolume = fVoxelizedDetectorZDivisionLog;
  fIsBuilt = true;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::UpdateVoxelizedGeometry()
{
  // Nothing happens if the voxelized geometry is not built. But parameters are properly set.
  if (!fIsBuilt) {
    return;
  }

  CalculateVoxelNumber();

  // Volume that will be deleted in order to update
  G4VPhysicalVolume* myVol;

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();

  myVol = store->GetVolume("VoxelizedDetectorXDivisionPhys");
  store->DeRegister(myVol);
  myVol = store->GetVolume("VoxelizedDetectorYDivisionPhys");
  store->DeRegister(myVol);
  myVol = store->GetVolume("VoxelizedDetectorZDivisionPhys");
  store->DeRegister(myVol);
  fVoxelizedDetectorXDivisionPhys =
    new G4PVReplica("VoxelizedDetectorXDivisionPhys", fVoxelizedDetectorXDivisionLog, fWorldLogical,
                    kXAxis, fVoxelNumberAlongX, fVoxelWidthX);

  fVoxelizedDetectorYDivisionPhys =
    new G4PVReplica("VoxelizedDetectorYDivisionPhys", fVoxelizedDetectorYDivisionLog,
                    fVoxelizedDetectorXDivisionPhys, kYAxis, fVoxelNumberAlongY, fVoxelWidthY);

  fVoxelizedDetectorZDivisionPhys =
    new G4PVReplica("VoxelizedDetectorZDivisionPhys", fVoxelizedDetectorZDivisionLog,
                    fVoxelizedDetectorYDivisionPhys, kZAxis, fVoxelNumberAlongZ, fVoxelWidthZ);

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructSD()
{
  G4String sensitiveDetectorName = "VoxelizedDetector";
  G4String HCname = "LETdata";

  SD* detectorSD = new SD(sensitiveDetectorName, HCname);
  G4SDManager::GetSDMpointer()->AddNewDetector(detectorSD);
  fSensitiveLogicalVolume->SetSensitiveDetector(detectorSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::Construct()
{
  ConstructVoxelizedDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::InitializeWorldPtr(G4VPhysicalVolume* pWorld)
{
  if (pWorld == nullptr)
    G4Exception("VoxelizedSensitiveDetector::InitializeWorldPtr", "WorldinitNull", FatalException,
                "Initializing Voxelization Class with a Null Pointer to World!");
  fWorldLogical = pWorld->GetLogicalVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Returns absolute voxel index given matrix indexes
G4int VoxelizedSensitiveDetector::GetThisVoxelNumber(G4int x, G4int y, G4int z) const
{
  G4int nz = GetVoxelNumberAlongZ();
  G4int ny = GetVoxelNumberAlongY();

  return z + nz * (y + ny * (x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
