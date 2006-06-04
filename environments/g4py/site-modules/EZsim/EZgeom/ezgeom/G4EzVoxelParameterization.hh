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
// $Id: G4EzVoxelParameterization.hh,v 1.2 2006-06-04 21:36:34 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   G4EzVoxelParameterization.hh
//
//                                         2005 Q
// ====================================================================
#ifndef G4EZ_VOXEL_PARAMETERIZATION_H
#define G4EZ_VOXEL_PARAMETERIZATION_H

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class G4EzVoxelParameterization : public G4VPVParameterisation {
protected:
  G4ThreeVector voxelSize;
  G4ThreeVector motherSize;
  G4int nvoxels[3];

public:
  G4EzVoxelParameterization(G4double ddx, G4double ddy, G4double ddz, 
			    G4int nx, G4int ny, G4int nz);

  ~G4EzVoxelParameterization();

  virtual void ComputeTransformation(const G4int copyNo,
				     G4VPhysicalVolume* physVol) const;
  
  virtual void ComputeDimensions(G4Box& aBox, const G4int copyNo,
                                 const G4VPhysicalVolume* physVol) const;
};

#endif
