// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorConstruction.cc,v 1.7 2000-01-17 10:30:17 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "MyDetectorConstruction.hh"

#include "MyDetectorMessenger.hh"
#include "MyCalorimeterSD.hh"
#include "MyTrackerSD.hh"
#include "MyCalorimeterHit.hh"
#include "MyTrackerHit.hh"
#include "MyCalorimeterHitsCollection.hh"
#include "MyTrackerHitsCollection.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
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
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Air,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);

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
    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(0.*cm,i*calPos,0.*cm)),
                      "calo_phys",calorimeter_log,experimentalHall_phys,
                      false,i);
  }

  //------------------------------ tracker tube
  G4Tubs * tracker_tube
    = new G4Tubs("tracker_tube",0.*cm,trackerRadius,trackerHight,
                 0.*deg,360.*deg);
  G4LogicalVolume * tracker_log
    = new G4LogicalVolume(tracker_tube,Ar,"tracker_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.*cm,trackerPos,0.*cm),
                    "tracker_phys",tracker_log,experimentalHall_phys,
                    false,0);

  //------------------------------ displaced box
  G4Box * undisplaced_box = new G4Box("undisplaced_box",30.*cm,50.*cm,70.*cm);
  G4DisplacedSolid* displaced_box = new G4DisplacedSolid
    ("displaced_box",undisplaced_box,
     G4Transform3D(G4RotationMatrix().rotateZ(20.*deg),
		   G4ThreeVector(200.*cm,0.,0.)));
  G4LogicalVolume * displaced_box_log
    = new G4LogicalVolume(displaced_box,Ar,"displaced_box_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.*cm,-200.*cm,0.*cm),
                    "displaced_box_phys",displaced_box_log,
		    experimentalHall_phys,
                    false,0);

  //-------------------------------------------- Boolean solids

  G4Tubs* cylinder1 = new G4Tubs("Cylinder #1",20*cm,50*cm,30*cm,0,2*M_PI);
  G4Box* box1 = new G4Box("Box #1",20*cm,30*cm,40*cm);
  G4Box* box2 = new G4Box("Box #2",10*cm,20*cm,30*cm);
  G4RotationMatrix* rm1 = new G4RotationMatrix;
  rm1->rotateZ(20*deg);
  G4RotationMatrix* rm2 = new G4RotationMatrix;
  rm2->rotateZ(60*deg);

  G4UnionSolid* cyl1Ubox1 =
    new G4UnionSolid("cylinder1-union-box1", cylinder1, box1,
		     rm1,G4ThreeVector(30.*cm,0.,0.));
  G4UnionSolid* cyl1Ubox1Ubox2 =
    new G4UnionSolid("cylinder1-union-box1-union-box2", cyl1Ubox1, box2,
		     rm2,G4ThreeVector(0.,40*cm,0.));
  G4LogicalVolume * union_log
    = new G4LogicalVolume(cyl1Ubox1Ubox2,Ar,"union_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(200.*cm,-50*cm,0.*cm),
                    "union_phys",union_log,experimentalHall_phys,
                    false,0);

  //----------- Tubes, replicas(!?) and daughter boxes

  const G4double eps = 10 * mm;
  const G4double alp = 10 * mrad;
  G4int iCopy;

  G4Tubs* tube = new G4Tubs("tube",20*cm,50*cm,30*cm,0,2*M_PI);
  G4LogicalVolume * tube_log
    = new G4LogicalVolume(tube,Ar,"tube_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(-200.*cm,0.,0.*cm),
                    "tube_phys",tube_log,experimentalHall_phys,
                    false,0);

  G4Tubs* divided_tube = new G4Tubs
    ("divided_tube",
     20*cm,50*cm,30*cm,0,M_PI/3.);
  G4LogicalVolume * divided_tube_log
    = new G4LogicalVolume(divided_tube,Ar,"divided_tube_L",0,0,0);
  /************ 
  new G4PVReplica("divided_tube_phys",divided_tube_log,tube_log,
  		  kPhi,6,M_PI/3.);
  *********/
  for (iCopy = 0; iCopy < 6; iCopy++) {
    new G4PVPlacement
      (G4Transform3D(G4RotationMatrix().rotateZ(iCopy*M_PI/3.),
		     G4ThreeVector()),
       divided_tube_log,"divided_tube_phys",tube_log,
       false,iCopy);
  }

  G4Tubs* divided_tube_inset = new G4Tubs
    ("divided_tube_inset",
     20*cm+2.*eps,50*cm-2.*eps,30*cm-2.*eps,2.*alp,M_PI/3.-4.*alp);
  G4LogicalVolume * divided_tube_inset_log
    = new G4LogicalVolume(divided_tube_inset,Ar,"divided_tube_inset_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(),
                    divided_tube_inset_log,"divided_tube_inset_phys",
		    divided_tube_log,false,0);

  G4Tubs* sub_divided_tube = new G4Tubs
    ("sub_divided_tube",
     20*cm+2.*eps,50*cm-2.*eps,30*cm-2.*eps,0.,M_PI/12.-alp);
  G4LogicalVolume * sub_divided_tube_log
    = new G4LogicalVolume(sub_divided_tube,Ar,"sub_divided_tube_L",0,0,0);
  /************ 
  new G4PVReplica("sub_divided_tube_phys",
  		  sub_divided_tube_log,divided_tube_inset_log,
		  kPhi,4,M_PI/12.-alp,2.*alp);
  ************/
  for (iCopy = 0; iCopy < 4; iCopy++) {
    new G4PVPlacement
      (G4Transform3D(G4RotationMatrix().rotateZ(iCopy*(M_PI/12.-alp)+2.*alp),
		     G4ThreeVector()),
       sub_divided_tube_log,"sub_divided_tube_phys",divided_tube_inset_log,
       false,iCopy);
  }
  
  G4Box* daughter_box = new G4Box("daughter_box",4.*cm,3.*cm,25.*cm);
  G4LogicalVolume * daughter_box_log
    = new G4LogicalVolume(daughter_box,Ar,"daughter_box_L",0,0,0);
  G4Box* grand_daughter_box = new G4Box("grand_daughter_box",1*cm,2*cm,5*cm);
  G4LogicalVolume * grand_daughter_box_log
    = new G4LogicalVolume(grand_daughter_box,Ar,"grand_daughter_box_L",0,0,0);
  G4Box* grand_daughter_box2 = new G4Box("grand_daughter_box2",1*cm,2*cm,5*cm);
  G4LogicalVolume * grand_daughter_box_log2
    = new G4LogicalVolume
    (grand_daughter_box2,Ar,"grand_daughter_box2_L",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(-2*cm,0.,0.),
                    grand_daughter_box_log,"grand_daughter_box_phys",
		    daughter_box_log,false,0);
  new G4PVPlacement(0,G4ThreeVector(2*cm,0.,0.),
                    grand_daughter_box_log2,"grand_daughter_box_phys2",
		    daughter_box_log,false,0);
  new G4PVPlacement(0,G4ThreeVector(40*cm,5.*cm,0.),
                    daughter_box_log,"daughter_box_phys",
		    sub_divided_tube_log,false,0);

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

  //-------------------------------------------- Sphere

  G4Sphere* PD_vol_crystal
      = new G4Sphere("Test",
                     100.*cm,           // inner radius
                     200.*cm,           // outer radius
                     0.,                // start phi
                     90.*degree,         // delta phi
                     0.,                // start theta
                     90.*degree          // delta theta
                     );


  G4LogicalVolume* PD_log_crystal
    = new G4LogicalVolume(PD_vol_crystal,Ar,"Test");

  G4RotationMatrix rm;

  G4PVPlacement* PD_physical
    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(200.*cm,200.*cm,0)),
                        "PD_physical", PD_log_crystal,
                        experimentalHall_phys,false,0);

  //-------------------------------------------- visualization attributes

  //  G4VisAttributes * experimentalHallVisAtt
  //      = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  //  experimentalHallVisAtt->SetForceWireframe(true);
  //  experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);
  experimentalHall_log -> SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes * calorimeterVisAtt
      = new G4VisAttributes(G4Colour(0.,0.,1.));
  //  calorimeterVisAtt->SetForceWireframe(true);
  calorimeter_log->SetVisAttributes(calorimeterVisAtt);

  G4VisAttributes * trackerVisAtt
    = new G4VisAttributes(G4Colour(0.,0.,1.));
  //  trackerVisAtt->SetForceWireframe(true);
  tracker_log->SetVisAttributes(trackerVisAtt);

  //------------------------------------------------------------------

  return experimentalHall_phys;
}
