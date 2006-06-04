//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4EzVoxelParameterization.cc,v 1.2 2006-06-04 21:36:34 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   G4EzVoxelParameterization.cc
//
//                                         2005 Q
// ====================================================================
#include "G4EzVoxelParameterization.hh"
#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////////////////////////////////////////
G4EzVoxelParameterization::G4EzVoxelParameterization
  (G4double ddx, G4double ddy, G4double ddz, G4int nx, G4int ny, G4int nz)
    : voxelSize(ddx, ddy, ddz), motherSize(ddx*nx, ddy*ny, ddz*nz)
//////////////////////////////////////////////////////////////////////////
{
  nvoxels[0]= nx;
  nvoxels[1]= ny;
  nvoxels[2]= nz;
}


///////////////////////////////////////////////////////
G4EzVoxelParameterization::~G4EzVoxelParameterization()
///////////////////////////////////////////////////////
{
}


///////////////////////////////////////////////////////////
void G4EzVoxelParameterization::ComputeTransformation
     (const G4int copyNo, G4VPhysicalVolume* physVol) const
///////////////////////////////////////////////////////////
{
  // copyNo-> index
  G4int nxy= nvoxels[0]*nvoxels[1];
  G4int iz= copyNo/nxy;
  G4int ixy= copyNo - iz*nxy;
  G4int iy= ixy/nvoxels[0];
  G4int ix= ixy-iy*nvoxels[0];

  // voxel position
  G4double px= voxelSize[0]*(ix+0.5) - motherSize[0]/2.;
  G4double py= voxelSize[1]*(iy+0.5) - motherSize[1]/2.;
  G4double pz= voxelSize[2]*(iz+0.5) - motherSize[2]/2.;
  
  physVol-> SetTranslation(G4ThreeVector(px, py, pz));

}


//////////////////////////////////////////////////////////////////////////////
void G4EzVoxelParameterization::ComputeDimensions
     (G4Box& aBox, const G4int copyNo, const G4VPhysicalVolume* physVol) const
//////////////////////////////////////////////////////////////////////////////
{
  return;
}

