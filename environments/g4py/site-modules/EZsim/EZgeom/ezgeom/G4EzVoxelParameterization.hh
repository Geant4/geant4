// $Id: G4EzVoxelParameterization.hh,v 1.1 2006-04-25 09:00:01 kmura Exp $
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
