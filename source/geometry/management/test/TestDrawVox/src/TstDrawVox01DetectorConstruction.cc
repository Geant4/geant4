
#include "TstDrawVox01DetectorConstruction.hh"

#include "TstDrawVox01DetectorMessenger.hh"

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
#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4TransportationManager.hh"


TstDrawVox01DetectorConstruction::TstDrawVox01DetectorConstruction()
:simpleBoxLog(NULL),simpleBoxDetector(NULL),honeycombDetector(NULL),
 detectorChoice(0),selectedMaterial(NULL),Air(NULL),Al(NULL),Pb(NULL)
{
  detectorMessenger = new TstDrawVox01DetectorMessenger(this);
  materialChoice = "Pb";
}

TstDrawVox01DetectorConstruction::~TstDrawVox01DetectorConstruction()
{
  delete detectorMessenger;
}

G4VPhysicalVolume* TstDrawVox01DetectorConstruction::Construct()
{
  if((!simpleBoxDetector)&&(!honeycombDetector))
  { ConstructDetectors(); }

  G4VPhysicalVolume* worldVol = NULL;
  switch(detectorChoice)
  { 
    case 1:
      worldVol = honeycombDetector; 
      break;
    default:
      worldVol = simpleBoxDetector;
  }
  return worldVol;
}

void TstDrawVox01DetectorConstruction::SwitchDetector()
{
  if((!simpleBoxDetector)&&(!honeycombDetector))
    ConstructDetectors();
  
  switch(detectorChoice)
  { 
    case 1:
      G4TransportationManager::GetTransportationManager()->
	GetNavigatorForTracking()->SetWorldVolume(honeycombDetector);
      break;
    default:
      G4TransportationManager::GetTransportationManager()->
	GetNavigatorForTracking()->SetWorldVolume(simpleBoxDetector);
  }
}

void TstDrawVox01DetectorConstruction::SelectDetector(G4String val)
{
  if(val=="Honeycomb") 
  { detectorChoice = 1; }
  else
  { detectorChoice = 0; }
  G4cout << "Now Detector is " << val << endl;
}

void TstDrawVox01DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << endl;
}

void TstDrawVox01DetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  if(!Air)
  {
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
    density = 1.29e-03*g/cm3;
    Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);
  }

  if(!Al)
  {
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    Al = new G4Material(name="Aluminium", z=13., a, density);
  }

  if(!Pb)
  {
    a = 207.19*g/mole;
    density = 11.35*g/cm3;
    Pb = new G4Material(name="Lead", z=82., a, density);
  }

  if(materialChoice=="Air")
  { selectedMaterial = Air; }
  else if(materialChoice=="Al")
  { selectedMaterial = Al; }
  else
  { selectedMaterial = Pb; }

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

void TstDrawVox01DetectorConstruction::ConstructDetectors()
{
  SelectMaterialPointer();

//--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------

//--------- simpleBoxDetector -------------------------------------

  G4Box * mySimpleBox = new G4Box("SBox",2000*cm, 2000*cm, 2000*cm);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

//--------- honeycombDetector -------------------------------------

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

  honeycombDetector = new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
  G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "CalPhys",
					       myCalLog,
                                               honeycombDetector,
					       false,0);

  G4String tName1("TPhys1");	// Allow all target physicals to share
				// same name (delayed copy)
  copyNo=0;
  for (j=1;j<=3;j++)
  {
    yTlate = -1000.0*cm - 40.0*cm + j*320.0*cm;
    for (i=1;i<= 4;i++)
    {
      xTlate = -1000.0*cm - 20.0*cm + i*160.0*cm - offset;
      G4PVPlacement *myTargetPhys
        =new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName1,myTargetLog,myCalPhys,false,copyNo++);
    }
  }
  for (j=1;j<=3;j++)
  {
    yTlate = -1000.0*cm - 80.0*cm + j*320.0*cm;
    for (i=1;i<= 4;i++)
    {
      xTlate = -1000.0*cm - 20.0*cm + i*160.0*cm;
      G4PVPlacement *myTargetPhys
        =new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName1,myTargetLog,myCalPhys,false,copyNo++);
    }
  }

//--------- Visualization attributes -------------------------------

  G4VisAttributes * simpleBoxVisAtt
    = new G4VisAttributes(G4Colour(.0,0.7,.0));
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxLog->SetVisAttributes(simpleBoxVisAtt);

  G4VisAttributes * experimentalHallVisAtt
    = new G4VisAttributes(G4Colour(0.3,0.,.0));
  experimentalHallVisAtt->SetVisibility(true);
  myWorldLog->SetVisAttributes(experimentalHallVisAtt);

  G4VisAttributes * calorimeterBoxVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  calorimeterBoxVisAtt->SetForceWireframe(true);
  calorimeterBoxVisAtt->SetVisibility(true);
  myCalLog->SetVisAttributes(calorimeterBoxVisAtt);

  G4VisAttributes * calorimeterTubeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.5,0.5));
  calorimeterTubeVisAtt->SetForceWireframe(true);
  calorimeterTubeVisAtt->SetVisibility(true);
  myTargetLog->SetVisAttributes(calorimeterTubeVisAtt);

//------------------------------------------------------------------

}

