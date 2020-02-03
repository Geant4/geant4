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
/// @file VoxelParam.cc
/// @brief Define voxel parameterization

#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "VoxelParam.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
VoxelParam::VoxelParam()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
VoxelParam::~VoxelParam()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void VoxelParam::ComputeTransformation(const G4int id,
                                       G4VPhysicalVolume* vol) const
{
  const G4int NX = 100;

  G4int iy = id / NX;
  G4int ix = id % NX;

  const G4double dxyz = 1.*mm;
  const G4double DXY = 10.*cm;

  G4double x0 = -DXY/2. + ix*dxyz;
  G4double y0 = -DXY/2. + iy*dxyz;

  vol-> SetTranslation(G4ThreeVector(x0,y0,0.));
  vol-> SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void VoxelParam::ComputeDimensions(G4Box& box, const G4int,
                                   const G4VPhysicalVolume*) const
{
  const G4double dxyz = 0.5*mm;

  box.SetXHalfLength(dxyz);
  box.SetYHalfLength(dxyz);
  box.SetZHalfLength(dxyz);
}
