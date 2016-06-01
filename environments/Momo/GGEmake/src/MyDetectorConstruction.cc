//------HeaderFile-
 #include "MyDetectorConstruction.hh"

#include "G4UnitsTable.hh"

#include "G4VUserDetectorConstruction.hh"

#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

MyDetectorConstruction::MyDetectorConstruction( )
{ ; }
MyDetectorConstruction::~MyDetectorConstruction( )
{ ; }
G4VPhysicalVolume* MyDetectorConstruction::Construct( )
{
// Elements
G4Element* elementN = new G4Element( "Nitrogen", "N", 7. , 14.00674*g/mole );
G4Element* elementO = new G4Element( "Oxygen", "O", 8. , 15.9994*g/mole );
G4Element* elementH = new G4Element( "Hydrogen", "H", 1. , 1.00794*g/mole );
G4Element* elementC = new G4Element( "Carbon", "C", 6. , 12.011*g/mole );

// Materials from Combination

G4Material* Air = new G4Material("Air",  1.205*g/cm3, 2, kStateGas, 293.0*kelvin, 1.0*atmosphere );
Air->AddElement( elementN, 3 );
Air->AddElement( elementO, 1 );
G4Material* Methane = new G4Material("Methane",  0.4241*g/cm3, 2, kStateLiquid, 111.7*kelvin, 1.0*atmosphere );
Methane->AddElement( elementH, 4 );
Methane->AddElement( elementC, 1 );

// Materials from Scratch

G4Material* Cu = new G4Material("Cu", 29, 63.546*g/mole, 8.96*g/cm3,kStateSolid, 293.0*kelvin, 1.0*pascal );
G4Material* Xeliquid = new G4Material("Xeliquid", 54, 131.29*g/mole, 3.52*g/cm3,kStateLiquid, 165.0*kelvin, 1.0*atmosphere );

// Visualization attributes


G4VisAttributes * blue= new G4VisAttributes( G4Colour(133,255,255 ));

G4VisAttributes * darkgray= new G4VisAttributes( G4Colour(75,84,87 ));

G4VisAttributes * trackred= new G4VisAttributes( G4Colour(234,84,87 ));

G4VisAttributes * calred= new G4VisAttributes( G4Colour(234,178,87 ));

// Logical  Volumes

G4Box *solidexpHall= new G4Box("solidexpHall", 5.0*m, 2.0*m, 2.0*m );
G4LogicalVolume * logicalexpHall = new G4LogicalVolume(solidexpHall, Air, "logicalexpHall" ,0,0,0);
logicalexpHall->SetVisAttributes(blue);
G4Box *solidcalblock= new G4Box("solidcalblock", 60.0*cm, 11.0*cm, 11.0*cm );
G4LogicalVolume * logicalcalblock = new G4LogicalVolume(solidcalblock, Xeliquid, "logicalcalblock" ,0,0,0);
logicalcalblock->SetVisAttributes(darkgray);
G4Tubs *solidtracker= new G4Tubs("solidtracker", 0.0*cm, 60.0*cm, 50.0*cm, 0.0*deg, 360.0*deg );
G4LogicalVolume * logicaltracker = new G4LogicalVolume(solidtracker, Methane, "logicaltracker" ,0,0,0);
logicaltracker->SetVisAttributes(trackred);
G4Box *solidcallayer= new G4Box("solidcallayer", 1.0*cm, 10.0*cm, 10.0*cm );
G4LogicalVolume * logicalcallayer = new G4LogicalVolume(solidcallayer, Cu, "logicalcallayer" ,0,0,0);
logicalcallayer->SetVisAttributes(calred);
G4Box *solidcalstack= new G4Box("solidcalstack", 60.0*cm, 11.0*cm, 68.0*cm );
G4LogicalVolume * logicalcalstack = new G4LogicalVolume(solidcalstack, Xeliquid, "logicalcalstack" ,0,0,0);
logicalcalstack->SetVisAttributes(darkgray);

// Physical Volumes -------Single Positioned Volume, Repeated Volumes, Replicas--------------------------- 


// Single Positioned Volumes 

G4VPhysicalVolume *  physicalexpHall= new G4PVPlacement(0,G4ThreeVector(0.0*m,0.0*mm,0.0*cm), "physicalexpHall",logicalexpHall,NULL,false, 0);

G4VPhysicalVolume *  physicaltracker= new G4PVPlacement(0,G4ThreeVector(-1.0*m,0.0*mm,0.0*cm), "physicaltracker",logicaltracker,physicalexpHall,false, 0);


// Repeated Volumes 

G4int copycallayerx;
copycallayerx=0;
for (G4int callayerx=1; callayerx<=20; callayerx++){
  G4double transcallayerx =-55.0*cm+5.0*cm*(callayerx-1);
  G4VPhysicalVolume * physicalcallayerx = new G4PVPlacement(0,G4ThreeVector(transcallayerx,0.0*cm,0.0*cm),  logicalcallayer,"physicalcallayerx",logicalcalblock,false,copycallayerx++);
}

G4int copystacky;
copystacky=0;
for (G4int stacky=1; stacky<=6; stacky++){
  G4double transstacky =-66.0*cm+22.0*cm*(stacky-1);
  G4VPhysicalVolume * physicalstacky = new G4PVPlacement(0,G4ThreeVector(0.0*cm, transstacky, 0.0*mm),  logicalcalblock,"physicalstacky",logicalcalstack,false,copystacky++);
}
G4int copystackz;
copystackz=0;
for (G4int stackz=1; stackz<=6; stackz++){
  G4double transstackz =-68.0*cm+22.0*cm*(stackz-1);
  G4VPhysicalVolume * physicalstackz = new G4PVPlacement(0,G4ThreeVector(80.0*cm,0.0*cm,transstackz), "physicalstackz", logicalcalstack, physicalexpHall,false,copystackz++);
}

// Replicas 


// return the physical World


 return physicalexpHall;
}
