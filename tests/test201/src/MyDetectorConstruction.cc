// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorConstruction.cc,v 1.4 2000-05-26 13:11:39 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "MyDetectorConstruction.hh"

#include "MyDetectorMessenger.hh"
#include "MyCalorimeterSD.hh"
#include "MyTrackerSD.hh"
#include "MyCalorimeterHit.hh"
#include "MyTrackerHit.hh"

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
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
  new MyDetectorMessenger(this);

  expHall_x = 600.*cm;
  expHall_y = 600.*cm;
  expHall_z = 600.*cm;

  calBox_x = 50.*cm;
  calBox_y = 50.*cm;
  calBox_z = 50.*cm;
  rotAngle = 30.*deg;
  calPos = 200.*cm;
  calMaterialName = "Pb";

  trackerRadius = 50.*cm;
  trackerHight = 100.*cm;
  trackerPos = -200.*cm;
}

MyDetectorConstruction::~MyDetectorConstruction()
{;}

G4VPhysicalVolume* MyDetectorConstruction::Construct()
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
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminum", z=13., a, density);

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);

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
  G4Material* calMat;
  if(calMaterialName=="Pb")
  { calMat = Pb; }
  else if(calMaterialName=="Al")
  { calMat = Al; }
  else if(calMaterialName=="Fe")
  { calMat = Fe; }
  else
  { calMat = Air; }
  G4LogicalVolume * calorimeter_log
    = new G4LogicalVolume(calorimeter_box,calMat,"calo_L",0,0,0);
  for(G4int i=0;i<3;i++)
  {
    G4RotationMatrix rm;
    rm.rotateZ(i*rotAngle);
    char s[64];
    sprintf(s,"calo_phys_%d",i);
    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(0.*cm,i*calPos,0.*cm)),
                      s,calorimeter_log,experimantalHall_phys,
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

  //------------------------------------------------ sensitive detectors

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String calorimeterSDname = "example2/calorimeter";
  MyCalorimeterSD * myCalorimeterSD = new MyCalorimeterSD( calorimeterSDname );
  SDman->AddNewDetector( myCalorimeterSD );
  calorimeter_log->SetSensitiveDetector( myCalorimeterSD );

  G4String trackerSDname = "example2/tracker";
  MyTrackerSD * myTrackerSD = new MyTrackerSD( trackerSDname );
  SDman->AddNewDetector( myTrackerSD );
  tracker_log->SetSensitiveDetector( myTrackerSD );

  //-------------------------------------------- visualization attributes

  //  G4VisAttributes * experimantalHallVisAtt
  //      = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  //  experimantalHallVisAtt->SetForceWireframe(true);
  //  experimantalHall_log->SetVisAttributes(experimantalHallVisAtt);
  experimantalHall_log -> SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes * calorimeterVisAtt
      = new G4VisAttributes(G4Colour(0.,0.,1.));
  //  calorimeterVisAtt->SetForceWireframe(true);
  calorimeter_log->SetVisAttributes(calorimeterVisAtt);

  G4VisAttributes * trackerVisAtt
    = new G4VisAttributes(G4Colour(0.,0.,1.));
  //  trackerVisAtt->SetForceWireframe(true);
  tracker_log->SetVisAttributes(trackerVisAtt);

  //------------------------------------------------------------------

  return experimantalHall_phys;
}

