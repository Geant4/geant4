
#include "MyDetectorConstruction.hh"

#include "MyTrackerSD.hh"

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

MyDetectorConstruction::MyDetectorConstruction()
{;}

MyDetectorConstruction::~MyDetectorConstruction()
{;}

G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
    G4cout << "  -- Tracker type sensitive detector" << endl;

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

  G4Box *myWorldBox= new G4Box("WBox",10000*cm, 10000*cm, 10000*cm);
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

  G4String tName1("TPhys1");
  copyNo=1;
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
  G4String tName2("TPhys2");
  copyNo=1;
  for (j=1;j<=26;j++)
  {
    yTlate = -1000.0*cm - 80.0*cm + j*80.0*cm;
    for (i=1;i<=50;i++)
    {
      xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm;
      G4PVPlacement *myTargetPhys
        =new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName2,myTargetLog,myCalPhys,false,copyNo++);
    }
  }

//--------- Sensitive detector -------------------------------------

 // A concrete object named "mySD", which corresponds to
 // the G4LogicalVolume object "myTargetLog"
    G4String sdName = "sample/tracker/mySD";
    MyTrackerSD* mySD = new MyTrackerSD( sdName );
    myTargetLog->SetSensitiveDetector( mySD );

 // "mySD" will be registrated
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector( mySD );

//------------------------------------------------------------------

  return myWorldPhys;
}

