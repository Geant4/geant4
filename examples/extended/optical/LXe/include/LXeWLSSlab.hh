#ifndef LXeWLSSlab_H
#define LXeWLSSlab_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "LXeDetectorConstruction.hh"

class LXeWLSSlab : public G4PVPlacement
{
public:
  LXeWLSSlab(G4RotationMatrix *pRot,
		const G4ThreeVector &tlate,
		G4LogicalVolume *pMotherLogical,
		G4bool pMany,
		G4int pCopyNo,
		LXeDetectorConstruction* c);
private:
  void CopyValues();
  
  LXeDetectorConstruction* constructor;

  G4bool updated;

  static G4LogicalVolume* ScintSlab_log;

  G4int nfibers;
  G4double scint_x;
  G4double scint_y;
  G4double scint_z;
  G4double slab_z;
};

#endif
