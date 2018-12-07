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
/// \file runAndEvent/RE02/src/RE02NestedPhantomParameterisation.cc
/// \brief Implementation of the RE02NestedPhantomParameterisation class
//
//
///////////////////////////////////////////////////////////////////////////////
#include "RE02NestedPhantomParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

//=======================================================================
// (RE02NestedPhantomParameterisation)
//
//  (Description)
//     Class for nested parameterisation.
//     This parameterisation handles material and transfomation of voxles.
//
//  T.Aso Created. Nov.2007.
//
////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02NestedPhantomParameterisation
::RE02NestedPhantomParameterisation(const G4ThreeVector& voxelSize,
                                    G4int nz,
                                    std::vector<G4Material*>& mat):
  G4VNestedParameterisation(),
  fdX(voxelSize.x()),fdY(voxelSize.y()),fdZ(voxelSize.z()),
  fNz(nz),fMat(mat)
{
  // Position of voxels. 
  // x and y positions are already defined in DetectorConstruction 
  // by using replicated volume. Here only we need to define is z positions
  // of voxles.
  fpZ.clear();
  G4double zp;
  for ( G4int iz = 0; iz < fNz; iz++){
    zp = (-fNz+1+2*iz)*fdZ;
    fpZ.push_back(zp);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02NestedPhantomParameterisation::~RE02NestedPhantomParameterisation(){
  fpZ.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Material assignment to geometry.
//
G4Material* RE02NestedPhantomParameterisation
::ComputeMaterial(G4VPhysicalVolume* /*currentVol*/, const G4int copyNo, 
                  const G4VTouchable* parentTouch)
{
  if(parentTouch==0) return fMat[0]; // protection for initialization and
                                     // vis at idle state
  // Copy number of voxels. 
  // Copy number of X and Y are obtained from replication number.
  // Copy nymber of Z is the copy number of current voxel.
  G4int ix = parentTouch->GetReplicaNumber(0);
  G4int iy = parentTouch->GetReplicaNumber(1);
  G4int iz = copyNo;
  // For demonstration purpose,a couple of materials are chosen alternately.
  G4Material* mat=0;
  if ( ix%2 == 0 && iy%2 == 0 && iz%2 == 0 ) mat = fMat[0];
  else mat = fMat[1];

  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Number of Materials
//  Material scanner is required for preparing physics tables and so on before 
//  stating simulation, so that G4 has to know number of materials.
G4int RE02NestedPhantomParameterisation::GetNumberOfMaterials() const{
  return fMat.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// GetMaterial
//  This is needed for material scanner and realizing geometry.
//
G4Material* RE02NestedPhantomParameterisation::GetMaterial(G4int i) const{
  return fMat[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Transformation of voxels.
//
void RE02NestedPhantomParameterisation
::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const{
  G4ThreeVector position(0.,0.,fpZ[copyNo]);
  physVol->SetTranslation(position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Dimensions are always same in this RE02 example.
//
void RE02NestedPhantomParameterisation
::ComputeDimensions(G4Box& box, const G4int, const G4VPhysicalVolume* ) const{
  box.SetXHalfLength(fdX);
  box.SetYHalfLength(fdY);
  box.SetZHalfLength(fdZ);
}
