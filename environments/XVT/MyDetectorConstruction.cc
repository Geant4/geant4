// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorConstruction.cc,v 1.2 1999-12-15 14:48:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyDetectorConstruction.hh"

#include "MyCalorimeterSD.cc"
#include "MyCalorimeterHit.cc"
#include "MyCalorimeterHitsCollection.cc"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

MyDetectorConstruction::MyDetectorConstruction()
{;}

MyDetectorConstruction::~MyDetectorConstruction()
{;}

G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
  G4cout << "Calorimeter volume Performance Test" << G4endl;

//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

//--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------

  G4double offset=22.5*cm, xTlate, yTlate;
  G4int i,j,copyNo;

  G4Box *myWorldBox= new G4Box("WBox",2000*cm, 2000*cm, 2000*cm);
  G4Box *myCalBox = new G4Box("CBox",1500*cm, 1500*cm, 1000*cm);
  G4Tubs *myTargetTube 
     = new G4Tubs("TTube",0*cm, 22.5*cm, 1000*cm, 0.*deg, 360.*deg);

  G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,Air,
						    "WLog", 0, 0, 0);
  G4LogicalVolume *myCalLog=new G4LogicalVolume(myCalBox,Al,
						  "CLog", 0, 0, 0);
  G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetTube,Pb,
						     "TLog", 0, 0, 0);

  G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
  G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "CalPhys",
					       myCalLog,
					       myWorldPhys,false,0);

  G4String tName1("TPhys1");	// Allow all target physicals to share
				// same name (delayed copy)
  copyNo=0;
  for (j=1;j<=25;j++)
  {
    yTlate = -1000.0*cm - 40.0*cm + j*80.0*cm;
    for (i=1;i<=50;i++)
    {
      xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm - offset;
      G4PVPlacement *myTargetPhys
        =new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName1,myTargetLog,myCalPhys,false,copyNo++);
    }
  }
  for (j=1;j<=26;j++)
  {
    yTlate = -1000.0*cm - 80.0*cm + j*80.0*cm;
    for (i=1;i<=50;i++)
    {
      xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm;
      G4PVPlacement *myTargetPhys
        =new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName1,myTargetLog,myCalPhys,false,copyNo++);
    }
  }

//--------- Sensitive detector -------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String calorimeterSDname = "example4/calorimeter";
  MyCalorimeterSD * myCalorimeterSD = new MyCalorimeterSD( calorimeterSDname );
  SDman->AddNewDetector( myCalorimeterSD );
  myTargetLog->SetSensitiveDetector(myCalorimeterSD);

//--------- Visualization attributes -------------------------------

  G4VisAttributes * experimantalHallVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  experimantalHallVisAtt->SetVisibility(false);
  myWorldLog->SetVisAttributes(experimantalHallVisAtt);

  G4VisAttributes * calorimeterBoxVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,0.5));
  calorimeterBoxVisAtt->SetForceWireframe(true);
  myCalLog->SetVisAttributes(calorimeterBoxVisAtt);

  G4VisAttributes * calorimeterTubeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.5,0.5));
  //  calorimeterTubeVisAtt->SetForceWireframe(true);
  myTargetLog->SetVisAttributes(calorimeterTubeVisAtt);

//------------------------------------------------------------------

  return myWorldPhys;
}

