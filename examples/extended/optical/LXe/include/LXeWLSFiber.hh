#ifndef LXeMainVolume_H
#define LXeMainVolume_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "LXeDetectorConstruction.hh"

class LXeWLSFiber : public G4PVPlacement
{
public:
  LXeWLSFiber(G4RotationMatrix *pRot,
		const G4ThreeVector &tlate,
		G4LogicalVolume *pMotherLogical,
		G4bool pMany,
		G4int pCopyNo,
		LXeDetectorConstruction* c);
private:

  void CopyValues();

  static G4LogicalVolume* clad2_log;

  G4bool updated; //does the fiber need to be rebuilt
  
  G4double fiber_rmin;    
  G4double fiber_rmax;    
  G4double fiber_z;
  G4double fiber_sphi;
  G4double fiber_ephi;

  G4double clad1_rmin;
  G4double clad1_rmax;    
  G4double clad1_z;
  G4double clad1_sphi;
  G4double clad1_ephi; 
  
  G4double clad2_rmin;
  G4double clad2_rmax;    
  G4double clad2_z;
  G4double clad2_sphi;
  G4double clad2_ephi;

  LXeDetectorConstruction* constructor;
};

#endif
