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
// $Id: DicomNestedPhantomParameterisation.cc,v 1.6 2009-11-10 18:48:46 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
#include "DicomNestedPhantomParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

DicomNestedPhantomParameterisation::
DicomNestedPhantomParameterisation(const G4ThreeVector& voxelSize,
                                         std::vector<G4Material*>& mat):
  G4VNestedParameterisation(), fdX(voxelSize.x()),
  fdY(voxelSize.y()), fdZ(voxelSize.z()), fMaterials(mat)
{
  // Position of voxels. 
  // x and y positions are already defined in DetectorConstruction by using
  // replicated volume. Here only we need to define is z positions of voxels.
}

DicomNestedPhantomParameterisation::~DicomNestedPhantomParameterisation()
{
}

void DicomNestedPhantomParameterisation::
SetNoVoxel( unsigned int nx, unsigned int ny, unsigned int nz )
{
  fnX = nx;
  fnY = ny;
  fnZ = nz;
}

//
// Material assignment to geometry.
//
G4Material* DicomNestedPhantomParameterisation::
ComputeMaterial(G4VPhysicalVolume*, const G4int copyNoZ, 
                                    const G4VTouchable* parentTouch)
{
  // protection for initialization and vis at idle state
  //
  if(parentTouch==0) return fMaterials[0];

  // Copy number of voxels. 
  // Copy number of X and Y are obtained from replication number.
  // Copy nymber of Z is the copy number of current voxel.
  G4int ix = parentTouch->GetReplicaNumber(0);
  G4int iy = parentTouch->GetReplicaNumber(1);
  G4int iz = copyNoZ;

  G4int copyNo = ix + fnX*iy + fnX*fnY*iz;

  unsigned int matIndex = GetMaterialIndex(copyNo);

  return fMaterials[ matIndex ];
}

//------------------------------------------------------------------
unsigned int DicomNestedPhantomParameterisation::
GetMaterialIndex( unsigned int copyNo ) const
{
  return *(fMaterialIndices+copyNo);
}

//
// Number of Materials
// Material scanner is required for preparing physics tables and so on before 
// starting simulation, so that G4 has to know number of materials.
//
G4int DicomNestedPhantomParameterisation::GetNumberOfMaterials() const
{
  return fMaterials.size();
}

//
// GetMaterial
//  This is needed for material scanner and realizing geometry.
//
G4Material* DicomNestedPhantomParameterisation::GetMaterial(G4int i) const
{
  return fMaterials[i];
}

//
// Transformation of voxels.
//
void DicomNestedPhantomParameterisation::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector position(0.,0.,(2*copyNo+1)*fdZ - fdZ*fnZ);
  physVol->SetTranslation(position);
}

//
// Dimensions are always same in this RE02 example.
//
void DicomNestedPhantomParameterisation::
ComputeDimensions( G4Box& box, const G4int, const G4VPhysicalVolume* ) const
{
  box.SetXHalfLength(fdX);
  box.SetYHalfLength(fdY);
  box.SetZHalfLength(fdZ);
}
