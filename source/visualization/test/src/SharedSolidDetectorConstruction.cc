#include "SharedSolidDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

SharedSolidDetectorConstruction::SharedSolidDetectorConstruction() {}

SharedSolidDetectorConstruction::~SharedSolidDetectorConstruction() {}

G4VPhysicalVolume* SharedSolidDetectorConstruction::Construct()
{
  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",10.*m,10.*m,10.*m);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,0,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);
  experimentalHall_log -> SetVisAttributes (G4VisAttributes::Invisible);

  G4Box* solidVoxel = new G4Box("Voxel",
                                3.6*mm/2,
                                3.6*mm/2,
                                3.6*mm/2);

  G4LogicalVolume* logicLung = new G4LogicalVolume(solidVoxel,
                                                     0,
                                                     "Lung");
  G4LogicalVolume* logicHeart = new G4LogicalVolume(solidVoxel,
                                                     0,
                                                     "Heart");

  new G4PVPlacement(0,
		    G4ThreeVector(0.,0.,-100.*mm),
		    logicLung,
		    "Lung",
		    experimentalHall_log,
		    false,
		    0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.,0.,+100.*mm),
		    logicHeart,
		    "Heart",
		    experimentalHall_log,
		    false,
		    0);
  logicLung->SetVisAttributes(new
			      G4VisAttributes(G4Colour(0.0,0.0,1.0)));
  logicHeart->SetVisAttributes(new
			       G4VisAttributes(G4Colour(1.0,0.0,0.0)));

  return experimentalHall_phys;
}
