//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//------------------------------ beam line along z axis---------------------

#include "MedLinacDetectorMessenger.hh"
#include "MedLinacDetectorConstruction.hh"
#include "MedLinacPhantomROGeometry.hh"
#include "MedLinacPhantomSD.hh"

#include "G4CSGSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4GeometryManager.hh"
#include "G4BooleanSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4VSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

MedLinacDetectorConstruction::MedLinacDetectorConstruction(G4String SDName)
  : phantomSD(0),phantomROGeometry(0),
    experimentalHall_log(0),
    //target_log(0), 
    vacuumBlock_log(0),
    collim_log(0), tracker_log(0),
    CylMinusCone_log(0),Filter_log(0),
    Mirror_log(0),
    JawY1_log(0),JawY2_log(0),
    JawX1_log(0),JawX2_log(0),
    Phantom_log(0),
    experimentalHall_phys(0),
    //target_phys(0), 
    vacuumBlock_phys(0),
    CylMinusCone_phys(0), Filter_phys(0),
    Mirror_phys(0),
    JawY1_phys(0), JawY2_phys(0),
    JawX1_phys(0), JawX2_phys(0),
    Phantom_phys(0)
{

  numberOfVoxelsAlongX=300;
  numberOfVoxelsAlongY=300;
  numberOfVoxelsAlongZ=300;

  ComputeDimVoxel();

  sensitiveDetectorName = SDName;

 // default parameter values
  JawX1Pos_x = -10. *cm;
  JawX2Pos_x =  10. *cm;
  JawY1Pos_y = -10. *cm;
  JawY2Pos_y =  10. *cm;

  detectorMessenger = new MedLinacDetectorMessenger(this);
 
}


MedLinacDetectorConstruction* MedLinacDetectorConstruction::instance = 0;

MedLinacDetectorConstruction* MedLinacDetectorConstruction::GetInstance(G4String SDName)
{
  if (instance == 0)
    {
      instance = new MedLinacDetectorConstruction(SDName);
     
    }
  return instance;
}



MedLinacDetectorConstruction::~MedLinacDetectorConstruction()
{

  delete detectorMessenger;

  if (phantomROGeometry) delete phantomROGeometry;
}


G4VPhysicalVolume* MedLinacDetectorConstruction::Construct()
{

  //---------rotation matrix first collimator and filter--------

  G4RotationMatrix*  rotateMatrix=new G4RotationMatrix();
  rotateMatrix->rotateX(180.0*deg);

  //---------rotation Mirror--------

  G4RotationMatrix*  rotateMirror=new G4RotationMatrix();
  rotateMirror->rotateX(35.0*deg);

  //---------------colors----------

  G4Colour  white   (1.0, 1.0, 1.0);
  G4Colour  grey    (0.5, 0.5, 0.5);
  //G4Colour  lgrey   (.75, .75, .75);
  //G4Colour  red     (1.0, 0.0, 0.0);
  //G4Colour  blue    (0.0, 0.0, 1.0);
  //G4Colour  cyan    (0.0, 1.0, 1.0);
  G4Colour  magenta (1.0, 0.0, 1.0); 
  //G4Colour  yellow  (1.0, 1.0, 0.0);
  G4Colour  lblue   (0.20, .50, .85);



  //-------------- materials------------------------------

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4int ncomponents;
  G4int natoms;
  G4String name;
  G4String symbol;
  G4double fractionmass;
  G4double massOfMole;                                   
  G4double temperature;                                             
  G4double pressure;                                              


  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen" ,symbol="N", z=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen",symbol="O", z=8., a);

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen",symbol="H", z=1., a);

  a = 49.29*g/mole;
  density = 2.6989*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air",density, ncomponents=2);
  Air->AddElement(elN, fractionmass=70*perCent);
  Air->AddElement(elO, fractionmass=30*perCent);

  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water",density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);


  massOfMole = 1.008*g/mole;
  density = 1.e-25*g/cm3;
  temperature = 2.73*kelvin;
  pressure = 3.e-18*pascal; 
  G4Material* Vacuum = new G4Material("interGalactic", z=1.,massOfMole, 
				      density, kStateGas,temperature, pressure);  
  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)
  //------------------------------ beam line along z axis---------------------

  G4double expHall_x = 3.0*m;
  G4double expHall_y = 3.0*m;
  G4double expHall_z = 3.0*m;
  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Air,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                                      "expHall",experimentalHall_log,0,false,0);
//------------------------target------------------------

  //G4double targetDim_x = 0.5*cm;
  //G4double targetDim_y = 0.5*cm;
  //G4double targetDim_z = 0.04*cm;
  //G4Box* target_box = new G4Box("target_box",targetDim_x,targetDim_y,targetDim_z);
  //target_log = new G4LogicalVolume(target_box,Pb,"target_log",0,0,0);
  //G4double targetPos_x = 0.0*m;
  //G4double targetPos_y = 0.0*m;
  //G4double targetPos_z = 99.96*cm;
  //target_phys = new G4PVPlacement(0,
  //         G4ThreeVector(targetPos_x,targetPos_y,targetPos_z),
  //         "target",target_log,experimentalHall_phys,false,0);
 
  //------------------------vacuum box------------------------
  G4double vacudim_x = 5.*cm;
  G4double vacudim_y = 5.*cm;
  G4double vacudim_z = 3.8*cm;
  G4Box* vacuumBlock_box = new G4Box("vacuumBlock_box",vacudim_x,
				     vacudim_y,vacudim_z);
  vacuumBlock_log = new G4LogicalVolume(vacuumBlock_box,
					Vacuum,"vacuumBlock_log",0,0,0);
  G4double vacublockPos_x = 0.0*m;
  G4double vacublockPos_y = 0.0*m;
  G4double vacublockPos_z = 96.1*cm;
  vacuumBlock_phys = new G4PVPlacement(0,
             G4ThreeVector(vacublockPos_x,vacublockPos_y,vacublockPos_z),
             "vacuBlock",vacuumBlock_log,experimentalHall_phys,false,0);

  //-------------------- the first collimator----------------

  G4double  pRmin1 = 0.*cm;
  G4double  pRmax1 = 0.3*cm;
  G4double  pRmin2 = 0.*cm;
  G4double  pRmax2 = 1.5*cm;
  G4double  hightOfTheCone =3.7*cm;
  G4double  startAngleOfTheCone = 0.*deg;
  G4double  spanningAngleOfTheCone = 360.*deg;

  G4Cons* collim_cone = new G4Cons("collim_cone",pRmin1,pRmax1,pRmin2,
				   pRmax2,hightOfTheCone,startAngleOfTheCone,
				   spanningAngleOfTheCone);
  collim_log = new G4LogicalVolume(collim_cone,Air,"collim_log",0,0,0);
 

  G4double innerRadiusOfTheTube = 0.*cm;
  G4double outerRadiusOfTheTube = 3.*cm;
  G4double hightOfTheTube = 3.7*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;
  G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
                                    outerRadiusOfTheTube,hightOfTheTube,
				    startAngleOfTheTube,spanningAngleOfTheTube);
  tracker_log = new G4LogicalVolume(tracker_tube,Pb,"tracker_log",0,0,0);


  G4SubtractionSolid* CylMinusCone = new G4SubtractionSolid("Cyl-Cone",
  							tracker_tube,collim_cone);
  CylMinusCone_log = new G4LogicalVolume(CylMinusCone,Pb,"CylminusCone_log",0,0,0);
  G4double CminusCPos_x = 0.*cm;
  G4double CminusCPos_y = 0.*cm;
  G4double CminusCPos_z = 96.1*cm;
  CylMinusCone_phys = new G4PVPlacement(rotateMatrix,
					G4ThreeVector(CminusCPos_x,CminusCPos_y,CminusCPos_z),
					"CylMinusCone",CylMinusCone_log,experimentalHall_phys,false,0);
  
  
  delete collim_log;
  delete tracker_log;

  //--------------filter--conic-------------

  G4double  FilterRmin1 = 0.*cm;
  G4double  FilterRmax1 = 0.1*cm;
  G4double  FilterRmin2 = 0.*cm;
  G4double  FilterRmax2 = 3.0*cm;
  G4double  FilterHightOfTheCone =0.45*cm;
  G4Cons* Filter = new G4Cons("Filter",FilterRmin1,FilterRmax1,FilterRmin2,
  			   FilterRmax2,FilterHightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  Filter_log = new G4LogicalVolume(Filter,Pb,"Filter_log",0,0,0);
  G4double FilterPos_x = 0.0*m;
  G4double FilterPos_y = 0.0*m;
  G4double FilterPos_z = 87.5*cm;
  Filter_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(FilterPos_x,FilterPos_y,FilterPos_z),
             "Filter",Filter_log,experimentalHall_phys,false,0);


  //----------------Mirror---------------

  G4double MirrorDim_x = 4.*cm;
  G4double MirrorDim_y = 4.*cm;
  G4double MirrorDim_z = 0.05*cm;
  G4Box* Mirror = new G4Box("Mirror",MirrorDim_x,
                                MirrorDim_y,MirrorDim_z);
  Mirror_log = new G4LogicalVolume(Mirror,Al,"Mirror_log",0,0,0);
  G4double MirrorPos_x = 0.0*cm;
  G4double MirrorPos_y = 0.0*cm;
  G4double MirrorPos_z = 72.0*cm;
  Mirror_phys = new G4PVPlacement(rotateMirror,
           G4ThreeVector(MirrorPos_x,MirrorPos_y,MirrorPos_z),
           "Mirror",Mirror_log,experimentalHall_phys,false,0);


  //----------------Jaw Y1 collimator---------------

  G4double JawY1Dim_x = 4.1*cm;
  G4double JawY1Dim_y = 4.1*cm;
  G4double JawY1Dim_z = 4.1*cm;
  G4Box* JawY1 = new G4Box("JawY1",JawY1Dim_x,
                                  JawY1Dim_y,JawY1Dim_z);
  JawY1_log = new G4LogicalVolume(JawY1,
                                             Pb,"JawY1_log",0,0,0);
  G4double JawY1Pos_x = 0.0*cm;
  //G4double JawY1Pos_y = -8.*cm;
  G4double JawY1Pos_z = 68.2*cm;
  JawY1_phys = new G4PVPlacement(0,
             G4ThreeVector(JawY1Pos_x,JawY1Pos_y,JawY1Pos_z),
             "JawY1",JawY1_log,experimentalHall_phys,false,0);

  //----------------Jaw Y2 collimator---------------

  G4double JawY2Dim_x = 4.1*cm;
  G4double JawY2Dim_y = 4.1*cm;
  G4double JawY2Dim_z = 4.1*cm;
  G4Box* JawY2 = new G4Box("JawY2",JawY2Dim_x,
                                  JawY2Dim_y,JawY2Dim_z);
  JawY2_log = new G4LogicalVolume(JawY2,
                                             Pb,"JawY1_log",0,0,0);
  G4double JawY2Pos_x = 0.0*cm;
  //G4double JawY2Pos_y = 8.*cm;
  G4double JawY2Pos_z = 68.2*cm;
  JawY2_phys = new G4PVPlacement(0,
             G4ThreeVector(JawY2Pos_x,JawY2Pos_y,JawY2Pos_z),
             "JawY2",JawY2_log,experimentalHall_phys,false,0);


 //----------------Jaw X1 collimator---------------

  G4double JawX1Dim_x = 4.1*cm;
  G4double JawX1Dim_y = 4.1*cm;
  G4double JawX1Dim_z = 4.1*cm;
  G4Box* JawX1 = new G4Box("JawX1",JawX1Dim_x,
                                  JawX1Dim_y,JawX1Dim_z);
  JawX1_log = new G4LogicalVolume(JawX1,
                                             Pb,"JawX1_log",0,0,0);

  // G4double JawX1Pos_x = -8.*cm;
  G4double JawX1Pos_y = 0.0*cm;
  G4double JawX1Pos_z = 59.4*cm;
  JawX1_phys = new G4PVPlacement(0,
             G4ThreeVector(JawX1Pos_x,JawX1Pos_y,JawX1Pos_z),
             "JawX1",JawX1_log,experimentalHall_phys,false,0);

  //----------------Jaw X2 collimator---------------

  G4double JawX2Dim_x = 4.1*cm;
  G4double JawX2Dim_y = 4.1*cm;
  G4double JawX2Dim_z = 4.1*cm;
  G4Box* JawX2 = new G4Box("JawX2",JawX2Dim_x,
                                  JawX2Dim_y,JawX2Dim_z);
  JawX2_log = new G4LogicalVolume(JawX2,
                                             Pb,"JawX1_log",0,0,0);
  //G4double JawX2Pos_x = 8.*cm;
  G4double JawX2Pos_y = 0.0*cm;
  G4double JawX2Pos_z = 59.4*cm;
  JawX2_phys = new G4PVPlacement(0,
             G4ThreeVector(JawX2Pos_x,JawX2Pos_y,JawX2Pos_z),
             "JawX2",JawX2_log,experimentalHall_phys,false,0);

  //----------------Phantom---------
  phantomDim_x = 15.*cm;
  phantomDim_y = 15.*cm;
  phantomDim_z = 15.*cm;
  G4Box* Phantom = new G4Box("Phantom",phantomDim_x,
                                  phantomDim_y,phantomDim_z);
  Phantom_log = new G4LogicalVolume(Phantom,
                                             H2O,"Phantom_log",0,0,0);
  G4double PhantomPos_x = 0.0*m;
  G4double PhantomPos_y = 0.0*m;
  G4double PhantomPos_z = 0.0*cm;
  Phantom_phys = new G4PVPlacement(0,
             G4ThreeVector(PhantomPos_x,PhantomPos_y,PhantomPos_z),
             "Phantom_phys",Phantom_log,experimentalHall_phys,false,0);

  PrintParameters();   

//--------- Visualization attributes -------------------------------


 
   G4VisAttributes* simpleH2OVisAtt= new G4VisAttributes(lblue);
   simpleH2OVisAtt->SetVisibility(true);
   simpleH2OVisAtt->SetForceSolid(true);
   //Phantom_log->SetVisAttributes(simpleH2OVisAtt);

   G4VisAttributes* simpleLeadWVisAtt= new G4VisAttributes(magenta);
   simpleLeadWVisAtt->SetVisibility(true);
   simpleLeadWVisAtt->SetForceWireframe(true);
   collim_log->SetVisAttributes(simpleLeadWVisAtt);
   CylMinusCone_log->SetVisAttributes(simpleLeadWVisAtt);

   G4VisAttributes* simpleLeadSVisAtt= new G4VisAttributes(magenta);
   simpleLeadSVisAtt->SetVisibility(true);
   simpleLeadSVisAtt->SetForceSolid(true);
   //target_log->SetVisAttributes(simpleLeadSVisAtt);
   Filter_log->SetVisAttributes(simpleLeadSVisAtt);
   //JawY1_log->SetVisAttributes(simpleLeadSVisAtt);
   //JawY2_log->SetVisAttributes(simpleLeadSVisAtt);
   //JawX1_log->SetVisAttributes(simpleLeadSVisAtt);
   //JawX2_log->SetVisAttributes(simpleLeadSVisAtt);

   G4VisAttributes* simpleAlVisAtt= new G4VisAttributes(grey);
   simpleAlVisAtt->SetVisibility(true);
   simpleAlVisAtt->SetForceSolid(true);
   //Mirror_log->SetVisAttributes(simpleAlVisAtt);

   G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(white);
   simpleWorldVisAtt->SetVisibility(false);
   experimentalHall_log ->SetVisAttributes(simpleWorldVisAtt);

   //to not see these components:

   Phantom_log->SetVisAttributes(simpleWorldVisAtt);
   JawY1_log->SetVisAttributes(simpleWorldVisAtt);
   JawY2_log->SetVisAttributes(simpleWorldVisAtt);
   JawX1_log->SetVisAttributes(simpleWorldVisAtt);
   JawX2_log->SetVisAttributes(simpleWorldVisAtt);
   Mirror_log->SetVisAttributes(simpleWorldVisAtt);
   //CylMinusCone_log->SetVisAttributes(simpleWorldVisAtt);
   //target_log->SetVisAttributes(simpleWorldVisAtt);
   //Filter_log->SetVisAttributes(simpleWorldVisAtt);
   //vacuumBlock_log->SetVisAttributes(simpleWorldVisAtt);


   //----------------------------------------------------------
   ConstructSensitiveDetector();

   return experimentalHall_phys;
}

void MedLinacDetectorConstruction::PrintParameters()
{
  G4cout <<"jaws1 x position "<< JawX1Pos_x/cm << " cm "<<G4endl ;
  G4cout <<"jaws2 x position "<< JawX2Pos_x/cm << " cm "<<G4endl ; 
  G4cout <<"jaws1 y position "<< JawY1Pos_y/cm << " cm "<<G4endl ;
  G4cout <<"jaws2 y position "<< JawY2Pos_y/cm << " cm "<<G4endl ; 
    
}



void MedLinacDetectorConstruction::SetJawX1Pos_x (G4double val)
{
  JawX1Pos_x = val;
  G4cout <<"==========================DetectorConstruction JawX1 "<< JawX1Pos_x/cm<<"cm"<<G4endl;
}


void MedLinacDetectorConstruction::SetJawX2Pos_x (G4double val)
{
  JawX2Pos_x = val;
  G4cout <<"==========================DetectorConstruction JawX2  "<< JawX2Pos_x/cm<<"cm"<<G4endl;
}

void MedLinacDetectorConstruction::SetJawY1Pos_y (G4double val)
{
  JawY1Pos_y = val;
  G4cout <<"==========================DetectorConstruction JawY1 "<< JawY1Pos_y/cm<<"cm"<<G4endl;
}


void MedLinacDetectorConstruction::SetJawY2Pos_y (G4double val)
{
  JawY2Pos_y = val;
  G4cout <<"==========================DetectorConstruction JawY2  "<< JawY2Pos_y/cm<<"cm"<<G4endl;
}



#include "G4RunManager.hh"

void MedLinacDetectorConstruction::UpdateGeometry()
{
  G4cout <<"@@@@@@@@@@@@@@@@@@@@@UpdateGeometry "<< G4endl;

  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}



void  MedLinacDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  if(!phantomSD)
  {
 
    phantomSD = new MedLinacPhantomSD(sensitiveDetectorName,numberOfVoxelsAlongX,numberOfVoxelsAlongY,numberOfVoxelsAlongZ);
    G4String ROGeometryName = "PhantomROGeometry";
    phantomROGeometry = new MedLinacPhantomROGeometry(ROGeometryName, phantomDim_x,phantomDim_y,phantomDim_z,
						  numberOfVoxelsAlongX,numberOfVoxelsAlongY,numberOfVoxelsAlongZ);


    phantomROGeometry->BuildROGeometry();
    phantomSD->SetROgeometry(phantomROGeometry);
    pSDManager->AddNewDetector(phantomSD);
    Phantom_log->SetSensitiveDetector(phantomSD);
  }
}














