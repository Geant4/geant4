
#include "ExE01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

ExE01DetectorConstruction::ExE01DetectorConstruction()
{
  expHall_x = 600.*cm;
  expHall_y = 600.*cm;
  expHall_z = 600.*cm;

  calBox_x = 50.*cm;
  calBox_y = 50.*cm;
  calBox_z = 50.*cm;
  rotAngle = 30.*deg;
  calPos = 200.*cm;

  trackerRadius = 50.*cm;
  trackerHight = 100.*cm;
  trackerPos = -200.*cm;
}

ExE01DetectorConstruction::~ExE01DetectorConstruction()
{;}

G4VPhysicalVolume* ExE01DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Lead = new G4Material(name="Lead", z=82., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall
  G4Box * experimantalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimantalHall_log
    = new G4LogicalVolume(experimantalHall_box,Air,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimantalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimantalHall_log,0,false,0);

  //------------------------------ calorimeter boxes
  G4Box * calorimeter_box
    = new G4Box("calorimeter_b",calBox_x,calBox_y,calBox_z);
  G4LogicalVolume * calorimeter_log
    = new G4LogicalVolume(calorimeter_box,Lead,"calo_L",0,0,0);
  for(G4int i=0;i<3;i++)
  {
    G4RotationMatrix rm;
    rm.rotateZ(i*rotAngle);
    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(0.*cm,i*calPos,0.*cm)),
                      "calo_phys",calorimeter_log,experimantalHall_phys,
                      false,i);
  }

  //------------------------------ tracker tube
  G4Tubs * tracker_tube
    = new G4Tubs("tracker_tube",0.*cm,trackerRadius,trackerHight,
                 0.*deg,360.*deg);
  G4LogicalVolume * tracker_log
    = new G4LogicalVolume(tracker_tube,Ar,"tracker_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.*cm,trackerPos,0.*cm),
                    "tracker_phys",tracker_log,experimantalHall_phys,
                    false,0);

  //------------------------------------------------------------------

  return experimantalHall_phys;
}

