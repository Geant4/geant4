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
// $Id: G4EzVoxelParameterization.hh 66892 2013-01-17 10:57:59Z gunter $
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
