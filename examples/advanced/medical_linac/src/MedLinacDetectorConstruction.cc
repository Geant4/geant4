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
// $Id: MedLinacDetectorConstruction.cc,v 1.2 2004-04-02 17:48:03 mpiergen Exp $
//
// Code developed by: M. Piergentili
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
    experimentalHall_log(0),targetA_log(0), 
    targetB_log(0), vacuumBlock_log(0),
    collim_log(0), tracker_log(0),
    CylMinusCone_log(0),layer1_log(0),layer2_log(0),
    layer3_log(0),layer4_log(0),layer5_log(0),layer6_log(0),
    layer7_log(0),layer8_log(0),layer9_log(0),layer10_log(0),
    layer11_log(0),layer12_log(0),layer13_log(0),layer14_log(0),
    layer15_log(0),layer16_log(0),layer17_log(0),layer18_log(0),
    layer19_log(0),layer20A_log(0),cone20_log(0), layer20_log(0),
    layer21_log(0),Mirror_log(0),
    JawY1_log(0),JawY2_log(0),
    JawX1_log(0),JawX2_log(0),
    Phantom_log(0),
    experimentalHall_phys(0),targetA_phys(0), 
    targetB_phys(0), vacuumBlock_phys(0),
    CylMinusCone_phys(0), layer1_phys(0),
    layer2_phys(0),layer3_phys(0),layer4_phys(0),layer5_phys(0),
    layer6_phys(0),layer7_phys(0),layer8_phys(0),layer9_phys(0),
    layer10_phys(0),layer11_phys(0),layer12_phys(0),layer13_phys(0),
    layer14_phys(0),layer15_phys(0),layer16_phys(0),layer17_phys(0),
    layer18_phys(0),layer19_phys(0),layer20_phys(0),
    layer21_phys(0),Mirror_phys(0),
    JawY1_phys(0), JawY2_phys(0),
    JawX1_phys(0), JawX2_phys(0),
    Phantom_phys(0)
{
  //G4double phantomDim_x = 15.*cm;
  //G4double phantomDim_y = 15.*cm;
  //G4double phantomDim_z = 15.*cm;

   //PhantomDimensionX = PhantomDim_x;
   //PhantomDimensionY = PhantomDim_y;
   //PhantomDimensionZ = PhantomDim_z;

  numberOfVoxelsAlongX=300;
  numberOfVoxelsAlongY=300;
  numberOfVoxelsAlongZ=300;

  ComputeDimVoxel();

  sensitiveDetectorName = SDName;

 // default parameter values

  fieldX1 = -10.*cm;
  fieldX2 =  10.*cm;
  fieldY1 = -10.*cm;
  fieldY2 =  10.*cm;


  detectorMessenger = new MedLinacDetectorMessenger(this);
 

  G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction fieldX1 "<<fieldX1/cm<<"cm"<<G4endl;

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

  //if (phantomROGeometry) delete phantomROGeometry;
  //if (phantomSD) delete phantomSD;
}


G4VPhysicalVolume* MedLinacDetectorConstruction::Construct()
{

  return ConstructGeom();
}


G4VPhysicalVolume* MedLinacDetectorConstruction::ConstructGeom ()
{
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


  //a = 207.19*g/mole;
  //density = 11.35*g/cm3;
  //G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

  a = 14.00674*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen" ,symbol="N", z=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen",symbol="O", z=8., a);

  a = 1.00794*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen",symbol="H", z=1., a);

  a = 12.011*g/mole;
  G4Element* elC = new G4Element(name="Carbon",symbol="C", z=6. , a );

  //density = 2.702*g/cm3;
  //a = 26.981539*g/mole;
  //G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air",density, ncomponents=2);
  Air->AddElement(elN, fractionmass=70*perCent);
  Air->AddElement(elO, fractionmass=30*perCent);

  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water",density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

  //density = 19.3*g/cm3;
  density = 18.*g/cm3;
  a = 183.85*g/mole;
  G4Material* W = new G4Material(name="Tungsten"  , z=74., a, density);

  density = 1.848*g/cm3;
  a = 9.012182*g/mole;
  G4Material* Be = new G4Material(name="Beryllium"  , z=4., a, density);


  density = 1.39*g/cm3;
  G4Material* Mylar = new G4Material(name="Mylar",density , ncomponents=3, kStateUndefined, 273.15*kelvin, 1.0*atmosphere );
  Mylar->AddElement( elH, 4 );
  Mylar->AddElement( elC, 5 );
  Mylar->AddElement( elO, 2 );

  density = 1.42*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton",density , ncomponents=4);
  Kapton->AddElement( elH, 02.6362 *perCent);
  Kapton->AddElement( elC, 69.1133 *perCent);
  Kapton->AddElement( elN, 07.3270 *perCent);
  Kapton->AddElement( elO, 20.9235 *perCent);


  massOfMole = 1.008*g/mole;
  density = 1.e-25*g/cm3;
  temperature = 2.73*kelvin;
  pressure = 3.e-18*pascal; 
  G4Material* Vacuum = new G4Material("interGalactic", z=1.,massOfMole, 
				      density, kStateGas,temperature, pressure);  



  //------------------------------ beam line along z axis------------------

 //---------rotation matrix first collimator and filter--------

  G4RotationMatrix*  rotateMatrix=new G4RotationMatrix();
  rotateMatrix->rotateX(180.0*deg);

  //---------rotation Mirror--------

  G4RotationMatrix*  rotateMirror=new G4RotationMatrix();
  rotateMirror->rotateY(35.0*deg);

  //---------jaws position--------

  //----1Y------------------------
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1Y  fieldX1 "<<fieldX1/cm<<"cm"<<G4endl;
  G4double thetaY1 = atan(-fieldY1/(100.*cm));// in rad
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1Y thetaY1  "<<thetaY1<<"rad"<<G4endl;
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1Y thetaY1 deg  "<<thetaY1/deg<<"deg"<<G4endl;
  G4double JawY1Pos_y = -((3.9*cm+28.*cm)*sin(thetaY1)+10.*cm*cos(thetaY1));
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1Y JawY1Pos_y  "<<JawY1Pos_y/cm<<"cm"<<G4endl;
  G4double JawY1Pos_z =100.*cm-((3.9*cm+28.*cm)*cos(thetaY1))+10.*cm*sin(thetaY1);
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1Y JawY1Pos_z  "<<JawY1Pos_z/cm<<"cm"<<G4endl;

  //----2Y-------------------------
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1Y fieldY2 "<<fieldY2/cm<<"cm"<<G4endl;
  G4double thetaY2 = atan(fieldY2/(100.*cm)); // in rad
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2Y thetaY2 rad  "<<thetaY2/rad<<"rad"<<G4endl;
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2Y thetaY2 deg  "<<thetaY2/deg<<"deg"<<G4endl;
  G4double JawY2Pos_y = (3.9*cm+28.*cm)*sin(thetaY2)+10.*cm*cos(thetaY2);
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2Y JawY2Pos_y  "<<JawY2Pos_y/cm<<"cm"<<G4endl;
  G4double JawY2Pos_z = 100.*cm-(3.9*cm+28.*cm)*cos(thetaY2)+10.*cm*sin(thetaY2);
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2Y JawY2Pos_z  "<<JawY2Pos_z/cm<<"cm"<<G4endl;

  //----1X-------------------------
  G4double thetaX1 = atan(-fieldX1/(100.*cm));
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1X thetaX1 rad  "<<thetaX1/rad<<"rad"<<G4endl;
  G4double JawX1Pos_x = -((45.8*cm)*sin(thetaX1)+10.*cm*cos(thetaX1));
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1X JawX1Pos_x  "<<JawX1Pos_x/cm<<"cm"<<G4endl;
  G4double JawX1Pos_z =100.*cm-(45.8*cm)*cos(thetaX1)+10.*cm*sin(thetaX1);
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw1X JawX1Pos_z  "<<JawX1Pos_z/cm<<"cm"<<G4endl;

   //----2X-------------------------
  G4double thetaX2 = atan(fieldX2/(100.*cm));
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2X thetaX2 rad  "<<thetaX2/rad<<"rad"<<G4endl;
  G4double JawX2Pos_x = (45.8*cm)*sin(thetaX2)+10.*cm*cos(thetaX2);
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2X JawX2Pos_x  "<<JawX2Pos_x/cm<<"cm"<<G4endl;
  G4double JawX2Pos_z = 100.*cm-(45.8*cm)*cos(thetaX2)+10.*cm*sin(thetaX2);
  //G4cout <<"_+_+_+_+_+_+_+_+_+_DetectorConstruction jaw2X JawX2Pos_z  "<<JawX2Pos_z/cm<<"cm"<<G4endl;

  //---------rotation jaw1Y--------

  G4RotationMatrix* rotatejaw1Y=new G4RotationMatrix();
  rotatejaw1Y->rotateX(thetaY1);

//---------rotation jaw2Y--------

  G4RotationMatrix* rotatejaw2Y=new G4RotationMatrix();
  rotatejaw2Y->rotateX(-thetaY2);

  //---------rotation jaw1X--------

  G4RotationMatrix*  rotatejaw1X=new G4RotationMatrix();
  rotatejaw1X->rotateY(-thetaX1);

  //---------rotation jaw2X--------

  G4RotationMatrix*  rotatejaw2X=new G4RotationMatrix();
  rotatejaw2X->rotateY(thetaX2);

  //---------------colors----------

  G4Colour  white   (1.0, 1.0, 1.0);
  G4Colour  grey    (0.5, 0.5, 0.5);
  //G4Colour  lgrey   (.75, .75, .75);
  G4Colour  red     (1.0, 0.0, 0.0);
  //G4Colour  blue    (0.0, 0.0, 1.0);
  G4Colour  cyan    (0.0, 1.0, 1.0);
  G4Colour  magenta (1.0, 0.0, 1.0); 
  G4Colour  yellow  (1.0, 1.0, 0.0);
  G4Colour  lblue   (0.20, .50, .85);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)
  //------------------------------ beam line along z axis---------------------  

//------------------------------ experimental hall (world volume)
  G4double expHall_x = 3.0*m;
  G4double expHall_y = 3.0*m;
  G4double expHall_z = 3.0*m;
  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Air,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                                      "expHall",experimentalHall_log,0,false,0);

  //--------------------window upper-------------------

  G4double windowUpDim_x = 4.*cm;
  G4double windowUpDim_y = 4.*cm;
  G4double windowUpDim_z = 0.0125*cm;
  G4Box* windowUp_box = new G4Box("windowUp_box",windowUpDim_x,windowUpDim_y,windowUpDim_z);
  windowUp_log = new G4LogicalVolume(windowUp_box,Be,"windowUp_log",0,0,0);
  G4double windowUpPos_x = 0.0*m;
  G4double windowUpPos_y = 0.0*m;
  G4double windowUpPos_z = 108.5225*cm;
  windowUp_phys = new G4PVPlacement(0,
            G4ThreeVector(windowUpPos_x,windowUpPos_y,windowUpPos_z),
            "windowUp",windowUp_log,experimentalHall_phys,false,0);


//------------------------target 6MV------------------------

  G4double targetADim_x = 0.3576*cm;
  G4double targetADim_y = 0.3576*cm;
  G4double targetADim_z = 0.04445*cm;
  G4Box* targetA_box = new G4Box("targetA_box",targetADim_x,targetADim_y,targetADim_z);
  targetA_log = new G4LogicalVolume(targetA_box,W,"targetA_log",0,0,0);
  G4double targetAPos_x = 0.0*m;
  G4double targetAPos_y = 0.0*m;
  G4double targetAPos_z = 99.95555*cm;
  targetA_phys = new G4PVPlacement(0,
            G4ThreeVector(targetAPos_x,targetAPos_y,targetAPos_z),
            "targetA",targetA_log,experimentalHall_phys,false,0);


  G4double targetBDim_x = 0.3576*cm;
  G4double targetBDim_y = 0.3576*cm;
  G4double targetBDim_z = 0.07874*cm;
  G4Box* targetB_box = new G4Box("targetB_box",targetBDim_x,targetBDim_y,targetBDim_z);
  targetB_log = new G4LogicalVolume(targetB_box,Cu,"targetB_log",0,0,0);
  G4double targetBPos_x = 0.0*m;
  G4double targetBPos_y = 0.0*m;
  G4double targetBPos_z = 99.83236*cm;
  targetB_phys = new G4PVPlacement(0,
            G4ThreeVector(targetBPos_x,targetBPos_y,targetBPos_z),
            "targetB",targetB_log,experimentalHall_phys,false,0);

 
  //------------------------vacuum box------------------------

  G4double vacudim_x = 1.*cm;
  G4double vacudim_y = 1.*cm;
  G4double vacudim_z = 0.67681*cm;
  G4Box* vacuumBlock_box = new G4Box("vacuumBlock_box",vacudim_x,
				     vacudim_y,vacudim_z);
  vacuumBlock_log = new G4LogicalVolume(vacuumBlock_box,
					Vacuum,"vacuumBlock_log",0,0,0);
  G4double vacublockPos_x = 0.0*m;
  G4double vacublockPos_y = 0.0*m;
  G4double vacublockPos_z = 99.07681*cm;
  vacuumBlock_phys = new G4PVPlacement(0,
             G4ThreeVector(vacublockPos_x,vacublockPos_y,vacublockPos_z),
             "vacuBlock",vacuumBlock_log,experimentalHall_phys,false,0);


  //-------------------- the first collimator upper----------------


  G4double innerRadiusOfTheTubeEx = 1.0*cm;
  G4double outerRadiusOfTheTubeEx = 5.*cm;
  G4double hightOfTheTubeEx = 4.0*cm;
  G4double startAngleOfTheTubeEx = 0.*deg;
  G4double spanningAngleOfTheTubeEx = 360.*deg;
  G4Tubs* UpperCollimator = new G4Tubs("UpperCollimator",innerRadiusOfTheTubeEx,
                                    outerRadiusOfTheTubeEx,hightOfTheTubeEx,
				    startAngleOfTheTubeEx,spanningAngleOfTheTubeEx);
  UpperCollimator_log = new G4LogicalVolume(UpperCollimator,W,"UpperCollimator_log",0,0,0);

  G4double UpperCollimatorPosX = 0.*cm;
  G4double UpperCollimatorPosY = 0.*cm;
  G4double UpperCollimatorPosZ = 102.4*cm;

  UpperCollimator_phys = new G4PVPlacement(0,
					   G4ThreeVector(UpperCollimatorPosX,UpperCollimatorPosY,
							 UpperCollimatorPosZ),"UpperCollimator",
					   UpperCollimator_log,experimentalHall_phys,false,0);
  //-------------------- the first collimator lower----------------

  G4double  pRmin1 = 0.*cm;
  G4double  pRmax1 = 0.3576423*cm;
  G4double  pRmin2 = 0.*cm;
  G4double  pRmax2 = 1.7435066*cm;
  G4double  hightOfTheCone =3.1*cm;
  G4double  startAngleOfTheCone = 0.*deg;
  G4double  spanningAngleOfTheCone = 360.*deg;

  G4Cons* collim_cone = new G4Cons("collim_cone",pRmin1,pRmax1,pRmin2,
				   pRmax2,hightOfTheCone,startAngleOfTheCone,
				   spanningAngleOfTheCone);
  collim_log = new G4LogicalVolume(collim_cone,Vacuum,"collim_log",0,0,0);


  G4double innerRadiusOfTheTube = 0.*cm;
  G4double outerRadiusOfTheTube = 5.*cm;
  G4double hightOfTheTube = 3.1*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;
  G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
                                    outerRadiusOfTheTube,hightOfTheTube,
				    startAngleOfTheTube,spanningAngleOfTheTube);
  tracker_log = new G4LogicalVolume(tracker_tube,W,"tracker_log",0,0,0);


  G4SubtractionSolid* CylMinusCone = new G4SubtractionSolid("Cyl-Cone",
  							tracker_tube,collim_cone);
  CylMinusCone_log = new G4LogicalVolume(CylMinusCone,W,"CylminusCone_log",0,0,0);
  G4double CminusCPos_x = 0.*cm;
  G4double CminusCPos_y = 0.*cm;
  G4double CminusCPos_z = 95.3*cm;
  CylMinusCone_phys = new G4PVPlacement(rotateMatrix,
					G4ThreeVector(CminusCPos_x,CminusCPos_y,CminusCPos_z),
					"CylMinusCone",CylMinusCone_log,experimentalHall_phys,false,0);
  
  
  //delete collim_log;
  //delete tracker_log;
  //------------------window lower-------------------


  G4double windowLowDim_x = 4.*cm;
  G4double windowLowDim_y = 4.*cm;
  G4double windowLowDim_z = 0.0127*cm;
  G4Box* windowLow_box = new G4Box("windowLow_box",windowLowDim_x,windowLowDim_y,windowLowDim_z);

  windowLow_log = new G4LogicalVolume(windowLow_box,Be,"windowLow_log",0,0,0);

  G4double windowLowPos_x = 0.0*m;
  G4double windowLowPos_y = 0.0*m;
  G4double windowLowPos_z = 90.9873*cm;

  windowLow_phys = new G4PVPlacement(0,
            G4ThreeVector(windowLowPos_x,windowLowPos_y,windowLowPos_z),
            "windowLow",windowLow_log,experimentalHall_phys,false,0);

  //--------------Flattening Filter---------------
 G4double  layerRmin1 = 0.*cm;
  G4double  layerRmin2 = 0.*cm;
  G4double layerPos_x = 0.0*m;
  G4double layerPos_y = 0.0*m;

  //-----------layer1-----------------------------
  G4double  layer1Rmax1 = 0.00000001*cm;
  G4double  layer1Rmax2 = 0.025*cm;
  G4double  layer1HightOfTheCone = 0.006*cm;
  G4Cons* layer1 = new G4Cons("layer1",layerRmin1,layer1Rmax1,layerRmin2,
  			   layer1Rmax2,layer1HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer1_log = new G4LogicalVolume(layer1,Cu,"layer1_log",0,0,0);

  G4double layer1Pos_z = 88.301*cm;
  layer1_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer1Pos_z),
             "layer1",layer1_log,experimentalHall_phys,false,0);

  //-----------layer2-----------------------------
  G4double  layer2Rmax1 = 0.025*cm;
  G4double  layer2Rmax2 = 0.050*cm;
  G4double  layer2HightOfTheCone = 0.0065*cm;
  G4Cons* layer2 = new G4Cons("layer2",layerRmin1,layer2Rmax1,layerRmin2,
  			   layer2Rmax2,layer2HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer2_log = new G4LogicalVolume(layer2,Cu,"layer2_log",0,0,0);

  G4double layer2Pos_z = 88.2885*cm;
  layer2_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer2Pos_z),
             "layer2",layer2_log,experimentalHall_phys,false,0);
   //-----------layer3-----------------------------
  G4double  layer3Rmax1 = 0.050*cm;
  G4double  layer3Rmax2 = 0.075*cm;
  G4double  layer3HightOfTheCone = 0.0075*cm;
  G4Cons* layer3 = new G4Cons("layer3",layerRmin1,layer3Rmax1,layerRmin2,
  			   layer3Rmax2,layer3HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer3_log = new G4LogicalVolume(layer3,Cu,"layer3_log",0,0,0);

  G4double layer3Pos_z = 88.2745*cm;
  layer3_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer3Pos_z),
             "layer3",layer3_log,experimentalHall_phys,false,0);
 //-----------layer4-----------------------------
  G4double  layer4Rmax1 = 0.075*cm;
  G4double  layer4Rmax2 = 0.100*cm;
  G4double  layer4HightOfTheCone = 0.009*cm;
  G4Cons* layer4 = new G4Cons("layer4",layerRmin1,layer4Rmax1,layerRmin2,
  			   layer4Rmax2,layer4HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer4_log = new G4LogicalVolume(layer4,Cu,"layer4_log",0,0,0);

  G4double layer4Pos_z = 88.258*cm;
  layer4_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer4Pos_z),
             "layer4",layer4_log,experimentalHall_phys,false,0);

  //-----------layer5-----------------------------
  G4double  layer5Rmax1 = 0.100*cm;
  G4double  layer5Rmax2 = 0.150*cm;
  G4double  layer5HightOfTheCone = 0.022*cm;
  G4Cons* layer5 = new G4Cons("layer5",layerRmin1,layer5Rmax1,layerRmin2,
  			   layer5Rmax2,layer5HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer5_log = new G4LogicalVolume(layer5,Cu,"layer5_log",0,0,0);

  G4double layer5Pos_z = 88.227*cm;
  layer5_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer5Pos_z),
             "layer5",layer5_log,experimentalHall_phys,false,0);

 //-----------layer6-----------------------------
  G4double  layer6Rmax1 = 0.150*cm;
  G4double  layer6Rmax2 = 0.200*cm;
  G4double  layer6HightOfTheCone = 0.024*cm;
  G4Cons* layer6 = new G4Cons("layer6",layerRmin1,layer6Rmax1,layerRmin2,
  			   layer6Rmax2,layer6HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer6_log = new G4LogicalVolume(layer6,Cu,"layer6_log",0,0,0);

  G4double layer6Pos_z = 88.181*cm;
  layer6_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer6Pos_z),
             "layer6",layer6_log,experimentalHall_phys,false,0);

  //-----------layer7-----------------------------
  G4double  layer7Rmax1 = 0.200*cm;
  G4double  layer7Rmax2 = 0.250*cm;
  G4double  layer7HightOfTheCone = 0.023*cm;
  G4Cons* layer7 = new G4Cons("layer7",layerRmin1,layer7Rmax1,layerRmin2,
  			   layer7Rmax2,layer7HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer7_log = new G4LogicalVolume(layer7,Cu,"layer7_log",0,0,0);

  G4double layer7Pos_z = 88.134*cm;
  layer7_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer7Pos_z),
             "layer7",layer7_log,experimentalHall_phys,false,0);
 //-----------layer8-----------------------------
  G4double  layer8Rmax1 = 0.250*cm;
  G4double  layer8Rmax2 = 0.300*cm;
  G4double  layer8HightOfTheCone = 0.023*cm;
  G4Cons* layer8 = new G4Cons("layer8",layerRmin1,layer8Rmax1,layerRmin2,
  			   layer8Rmax2,layer8HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer8_log = new G4LogicalVolume(layer8,Cu,"layer8_log",0,0,0);

  G4double layer8Pos_z = 88.088*cm;
  layer8_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer8Pos_z),
             "layer8",layer8_log,experimentalHall_phys,false,0);

  //-----------layer9-----------------------------
  G4double  layer9Rmax1 = 0.300*cm;
  G4double  layer9Rmax2 = 0.350*cm;
  G4double  layer9HightOfTheCone = 0.0225*cm;
  G4Cons* layer9 = new G4Cons("layer9",layerRmin1,layer9Rmax1,layerRmin2,
  			   layer9Rmax2,layer9HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer9_log = new G4LogicalVolume(layer9,Cu,"layer9_log",0,0,0);

  G4double layer9Pos_z = 88.0425*cm;
  layer9_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer9Pos_z),
             "layer9",layer9_log,experimentalHall_phys,false,0);
  //-----------layer10-----------------------------
  G4double  layer10Rmax1 = 0.350*cm;
  G4double  layer10Rmax2 = 0.400*cm;
  G4double  layer10HightOfTheCone = 0.022*cm;
  G4Cons* layer10 = new G4Cons("layer10",layerRmin1,layer10Rmax1,layerRmin2,
  			   layer10Rmax2,layer10HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer10_log = new G4LogicalVolume(layer10,Cu,"layer10_log",0,0,0);

  G4double layer10Pos_z = 87.998*cm;
  layer10_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer10Pos_z),
             "layer10",layer10_log,experimentalHall_phys,false,0);

  //-----------layer11-----------------------------
  G4double  layer11Rmax1 = 0.400*cm;
  G4double  layer11Rmax2 = 0.500*cm;
  G4double  layer11HightOfTheCone = 0.0425*cm;
  G4Cons* layer11 = new G4Cons("layer11",layerRmin1,layer11Rmax1,layerRmin2,
  			   layer11Rmax2,layer11HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer11_log = new G4LogicalVolume(layer11,Cu,"layer11_log",0,0,0);

  G4double layer11Pos_z = 87.9335*cm;
  layer11_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer11Pos_z),
             "layer11",layer11_log,experimentalHall_phys,false,0);
  //-----------layer12-----------------------------
  G4double  layer12Rmax1 = 0.500*cm;
  G4double  layer12Rmax2 = 0.600*cm;
  G4double  layer12HightOfTheCone = 0.0395*cm;
  G4Cons* layer12 = new G4Cons("layer12",layerRmin1,layer12Rmax1,layerRmin2,
  			   layer12Rmax2,layer12HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer12_log = new G4LogicalVolume(layer12,Cu,"layer12_log",0,0,0);

  G4double layer12Pos_z = 87.8515*cm;
  layer12_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer12Pos_z),
             "layer12",layer12_log,experimentalHall_phys,false,0);

  //-----------layer13-----------------------------
  G4double  layer13Rmax1 = 0.600*cm;
  G4double  layer13Rmax2 = 0.700*cm;
  G4double  layer13HightOfTheCone = 0.038*cm;
  G4Cons* layer13 = new G4Cons("layer13",layerRmin1,layer13Rmax1,layerRmin2,
  			   layer13Rmax2,layer13HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer13_log = new G4LogicalVolume(layer13,Cu,"layer13_log",0,0,0);

  G4double layer13Pos_z = 87.774*cm;
  layer13_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer13Pos_z),
             "layer13",layer13_log,experimentalHall_phys,false,0);
 //-----------layer14-----------------------------
  G4double  layer14Rmax1 = 0.700*cm;
  G4double  layer14Rmax2 = 0.800*cm;
  G4double  layer14HightOfTheCone = 0.0335*cm;
  G4Cons* layer14 = new G4Cons("layer14",layerRmin1,layer14Rmax1,layerRmin2,
  			   layer14Rmax2,layer14HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer14_log = new G4LogicalVolume(layer14,Cu,"layer14_log",0,0,0);

  G4double layer14Pos_z = 87.7025*cm;
  layer14_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer14Pos_z),
             "layer14",layer14_log,experimentalHall_phys,false,0);

  //-----------layer15-----------------------------
  G4double  layer15Rmax1 = 0.800*cm;
  G4double  layer15Rmax2 = 0.900*cm;
  G4double  layer15HightOfTheCone = 0.0325*cm;
  G4Cons* layer15 = new G4Cons("layer15",layerRmin1,layer15Rmax1,layerRmin2,
  			   layer15Rmax2,layer15HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer15_log = new G4LogicalVolume(layer15,Cu,"layer15_log",0,0,0);

  G4double layer15Pos_z = 87.6365*cm;
  layer15_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer15Pos_z),
             "layer15",layer15_log,experimentalHall_phys,false,0);
  //-----------layer16-----------------------------
  G4double  layer16Rmax1 = 0.900*cm;
  G4double  layer16Rmax2 = 1.000*cm;
  G4double  layer16HightOfTheCone = 0.028*cm;
  G4Cons* layer16 = new G4Cons("layer16",layerRmin1,layer16Rmax1,layerRmin2,
  			   layer16Rmax2,layer16HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer16_log = new G4LogicalVolume(layer16,Cu,"layer16_log",0,0,0);

  G4double layer16Pos_z = 87.576*cm;
  layer16_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer16Pos_z),
             "layer16",layer16_log,experimentalHall_phys,false,0);

  //-----------layer17-----------------------------
  G4double  layer17Rmax1 = 1.000*cm;
  G4double  layer17Rmax2 = 1.100*cm;
  G4double  layer17HightOfTheCone = 0.0275*cm;
  G4Cons* layer17 = new G4Cons("layer17",layerRmin1,layer17Rmax1,layerRmin2,
  			   layer17Rmax2,layer17HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer17_log = new G4LogicalVolume(layer17,Cu,"layer17_log",0,0,0);

  G4double layer17Pos_z = 87.5205*cm;
  layer17_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer17Pos_z),
             "layer17",layer17_log,experimentalHall_phys,false,0);

  //-----------layer18-----------------------------
  G4double  layer18Rmax1 = 1.100*cm;
  G4double  layer18Rmax2 = 1.205*cm;
  G4double  layer18HightOfTheCone = 0.019*cm;
  G4Cons* layer18 = new G4Cons("layer18",layerRmin1,layer18Rmax1,layerRmin2,
  			   layer18Rmax2,layer18HightOfTheCone,startAngleOfTheCone,
  			   spanningAngleOfTheCone);

  layer18_log = new G4LogicalVolume(layer18,Cu,"layer18_log",0,0,0);

  G4double layer18Pos_z = 87.474*cm;
  layer18_phys = new G4PVPlacement(rotateMatrix,
             G4ThreeVector(layerPos_x,layerPos_y,layer18Pos_z),
             "layer18",layer18_log,experimentalHall_phys,false,0);

  //-----------layer19-----------------------------
  G4double innerRadiusOfTheLayer19 = 0.0*cm;
  G4double outerRadiusOfTheLayer19 = 1.325*cm;
  G4double hightOfTheLayer19 = 0.040*cm;
  G4double startAngleOfTheLayer19 = 0.*deg;
  G4double spanningAngleOfTheLayer19 = 360.*deg;
  G4Tubs* layer19 = new G4Tubs("layer19",innerRadiusOfTheLayer19,
                                    outerRadiusOfTheLayer19,hightOfTheLayer19,
				    startAngleOfTheLayer19,spanningAngleOfTheLayer19);
  layer19_log = new G4LogicalVolume(layer19,Cu,"layer19_log",0,0,0);

  G4double layer19PosX = 0.*cm;
  G4double layer19PosY = 0.*cm;
  G4double layer19PosZ = 87.415*cm;

  layer19_phys = new G4PVPlacement(0,
					   G4ThreeVector(layer19PosX,layer19PosY,
							 layer19PosZ),"layer19",
					   layer19_log,experimentalHall_phys,false,0);
  //-----------layer20-----------------------------
  G4double innerRadiusOfTheLayer20A = 1.300*cm;
  G4double outerRadiusOfTheLayer20A = 1.325*cm;
  G4double hightOfTheLayer20A = 0.0225*cm;
  G4double startAngleOfTheLayer20A = 0.*deg;
  G4double spanningAngleOfTheLayer20A = 360.*deg;
  G4Tubs* layer20A = new G4Tubs("layer20A",innerRadiusOfTheLayer20A,
                                    outerRadiusOfTheLayer20A,hightOfTheLayer20A,
				    startAngleOfTheLayer20A,spanningAngleOfTheLayer20A);
  layer20A_log = new G4LogicalVolume(layer20A,Cu,"layer20A_log",0,0,0);

  G4double  cone20Rmin1 = 0.*cm;
  G4double  cone20Rmax1 = 1.325*cm;
  G4double  cone20Rmin2 = 0.*cm;
  G4double  cone20Rmax2 = 1.300*cm;
  G4double  hightOfThecone20 =0.0225*cm;
  G4double  startAngleOfThecone20 = 0.*deg;
  G4double  spanningAngleOfThecone20 = 360.*deg;

  G4Cons* cone20 = new G4Cons("cone20",cone20Rmin1,cone20Rmax1,cone20Rmin2,
				   cone20Rmax2,hightOfThecone20,startAngleOfThecone20,
				   spanningAngleOfThecone20);
  cone20_log = new G4LogicalVolume(cone20,Air,"cone20_log",0,0,0);

  G4SubtractionSolid* layer20 = new G4SubtractionSolid("layer20",layer20A,cone20);
  layer20_log = new G4LogicalVolume(layer20,Cu,"layer20_log",0,0,0);

  G4double layer20PosX = 0.*cm;
  G4double layer20PosY = 0.*cm;
  G4double layer20PosZ = 87.4775*cm;
  layer20_phys = new G4PVPlacement(0,
					   G4ThreeVector(layer20PosX,layer20PosY,
							 layer20PosZ),"layer20",
				   layer20_log,experimentalHall_phys,false,0);

  //-----------layer21-----------------------------
  G4double innerRadiusOfTheLayer21 = 1.325*cm;
  G4double outerRadiusOfTheLayer21 = 1.500*cm;
  G4double hightOfTheLayer21 = 0.0625*cm;
  G4double startAngleOfTheLayer21 = 0.*deg;
  G4double spanningAngleOfTheLayer21 = 360.*deg;
  G4Tubs* layer21 = new G4Tubs("layer21",innerRadiusOfTheLayer21,
                                    outerRadiusOfTheLayer21,hightOfTheLayer21,
				    startAngleOfTheLayer21,spanningAngleOfTheLayer21);
  layer21_log = new G4LogicalVolume(layer21,Cu,"layer21_log",0,0,0);

  G4double layer21PosX = 0.*cm;
  G4double layer21PosY = 0.*cm;
  G4double layer21PosZ = 87.4375*cm;

  layer21_phys = new G4PVPlacement(0,
					   G4ThreeVector(layer21PosX,layer21PosY,
							 layer21PosZ),"layer19",
					   layer21_log,experimentalHall_phys,false,0);



  //---------------Ion chamber------------------------

  G4double innerRadiusOfTheIonChamber = 0.*cm;
  G4double outerRadiusOfTheIonChamber = 4.7625*cm;
  G4double hightOfTheWindows = 0.0127*cm;
  G4double hightOfTheSignalPlates = 0.00254*cm;
  G4double startAngleOfTheIonChamber = 0.*deg;
  G4double spanningAngleOfTheIonChamber = 360.*deg;
  G4double IonChamberPos_x = 0.0*cm;
  G4double IonChamberPos_y = 0.0*cm;


  G4Tubs* Window = new G4Tubs("Window",innerRadiusOfTheIonChamber,
                                    outerRadiusOfTheIonChamber,hightOfTheWindows,
			      startAngleOfTheIonChamber,spanningAngleOfTheIonChamber);

  Window_log = new G4LogicalVolume(Window,Kapton,"Window_log",0,0,0);


  G4Tubs* SignalPlate = new G4Tubs("SignalPlate",innerRadiusOfTheIonChamber,
                                    outerRadiusOfTheIonChamber,hightOfTheSignalPlates,
			      startAngleOfTheIonChamber,spanningAngleOfTheIonChamber);

  SignalPlate_log = new G4LogicalVolume(SignalPlate ,Kapton,"SignalPlate_log",0,0,0);

  //-----Ion chamber----window 1---------------------

  G4double Window1Pos_z = 85.165*cm;
  Window1_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,Window1Pos_z),
           "Window",Window_log,experimentalHall_phys,false,0);

 //------Ion chamber---Signal Plate 1---------------------

  G4double SignalPlate1Pos_z = 84.927*cm;
  SignalPlate1_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate1Pos_z),
           "SignalPlate",SignalPlate_log,experimentalHall_phys,false,0);

 //------Ion chamber---Signal Plate 2---------------------

  G4double SignalPlate2Pos_z = 84.688*cm;
  SignalPlate2_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate2Pos_z),
           "SignalPlate",SignalPlate_log,experimentalHall_phys,false,0);

 //------Ion chamber---window 2---------------------

  G4double Window2Pos_z = 84.45*cm;
  Window2_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,Window2Pos_z),
           "Window",Window_log,experimentalHall_phys,false,0);

 //------Ion chamber---Signal Plate 3---------------------

  G4double SignalPlate3Pos_z = 84.212*cm;
  SignalPlate3_phys = new G4PVPlacement(0,
					G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate3Pos_z),
					"SignalPlate",SignalPlate_log,experimentalHall_phys,false,0);

 //------Ion chamber---Signal Plate 4---------------------

  G4double SignalPlate4Pos_z = 83.973*cm;
  SignalPlate4_phys = new G4PVPlacement(0,
					G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate4Pos_z),
					"SignalPlate",SignalPlate_log,experimentalHall_phys,false,0);

 //------Ion chamber---window 3---------------------

  G4double Window3Pos_z = 83.735*cm;
  Window3_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,Window3Pos_z),
           "Window",Window_log,experimentalHall_phys,false,0);


  //----------------Mirror---------------------------

  G4double innerRadiusOfTheMirror = 0.*cm;
  G4double outerRadiusOfTheMirror = 6.*cm;
  G4double hightOfTheMirror = 0.00254*cm;
  G4double startAngleOfTheMirror = 0.*deg;
  G4double spanningAngleOfTheMirror = 360.*deg;
  G4Tubs* Mirror = new G4Tubs("Mirror",innerRadiusOfTheMirror,
                                    outerRadiusOfTheMirror,hightOfTheMirror,
			      startAngleOfTheMirror,spanningAngleOfTheMirror);

  Mirror_log = new G4LogicalVolume(Mirror,Mylar,"Mirror_log",0,0,0);
  G4double MirrorPos_x = 0.0*cm;
  G4double MirrorPos_y = 0.0*cm;
  G4double MirrorPos_z = 78.0*cm;
  Mirror_phys = new G4PVPlacement(rotateMirror,
           G4ThreeVector(MirrorPos_x,MirrorPos_y,MirrorPos_z),
           "Mirror",Mirror_log,experimentalHall_phys,false,0);


  //----------------Jaw Y1 collimator---------------

  G4double JawY1Dim_x = 25.0*cm;
  G4double JawY1Dim_y = 10.0*cm;
  G4double JawY1Dim_z = 3.9*cm;
  G4Box* JawY1 = new G4Box("JawY1",JawY1Dim_x,
                                  JawY1Dim_y,JawY1Dim_z);
  JawY1_log = new G4LogicalVolume(JawY1,
                                             W,"JawY1_log",0,0,0);



  G4double JawY1Pos_x = 0.0*cm;
  //G4double JawY1Pos_y = fieldY1;
  //G4double JawY1Pos_z = 68.1*cm;
  JawY1_phys = new G4PVPlacement(rotatejaw1Y,
             G4ThreeVector(JawY1Pos_x,JawY1Pos_y,JawY1Pos_z),
             "JawY1",JawY1_log,experimentalHall_phys,false,0);

  //----------------Jaw Y2 collimator---------------

  G4double JawY2Dim_x = 25.0*cm;
  G4double JawY2Dim_y = 10.0*cm;
  G4double JawY2Dim_z = 3.9*cm;
  G4Box* JawY2 = new G4Box("JawY2",JawY2Dim_x,
                                  JawY2Dim_y,JawY2Dim_z);
  JawY2_log = new G4LogicalVolume(JawY2,
                                             W,"JawY1_log",0,0,0);
  G4double JawY2Pos_x = 0.0*cm;
  //G4double JawY2Pos_y = fieldY2;
  //G4double JawY2Pos_z = 68.1*cm;
  JawY2_phys = new G4PVPlacement(rotatejaw2Y,
             G4ThreeVector(JawY2Pos_x,JawY2Pos_y,JawY2Pos_z),
             "JawY2",JawY2_log,experimentalHall_phys,false,0);


 //----------------Jaw X1 collimator---------------

  G4double JawX1Dim_x = 10.0*cm;
  G4double JawX1Dim_y = 25.0*cm;
  G4double JawX1Dim_z = 3.9*cm;
  G4Box* JawX1 = new G4Box("JawX1",JawX1Dim_x,
                                 JawX1Dim_y,JawX1Dim_z);
  JawX1_log = new G4LogicalVolume(JawX1,
                                             W,"JawX1_log",0,0,0);
  G4cout <<"_________________DetectorConstruction fieldX1 "<< fieldX1/cm << G4endl;

  //G4double JawX1Pos_x = fieldX1;
  G4double JawX1Pos_y = 0.0*cm;
  //G4double JawX1Pos_z = 54.2*cm;
  JawX1_phys = new G4PVPlacement(rotatejaw1X,
             G4ThreeVector(JawX1Pos_x,JawX1Pos_y,JawX1Pos_z),
             "JawX1",JawX1_log,experimentalHall_phys,false,0);

  //----------------Jaw X2 collimator---------------

  G4double JawX2Dim_x = 10.0*cm;
  G4double JawX2Dim_y = 25.0*cm;
  G4double JawX2Dim_z = 3.9*cm;
  G4Box* JawX2 = new G4Box("JawX2",JawX2Dim_x,
                                  JawX2Dim_y,JawX2Dim_z);
  JawX2_log = new G4LogicalVolume(JawX2,
                                             W,"JawX1_log",0,0,0);
  //G4double JawX2Pos_x = fieldX2;
  G4double JawX2Pos_y = 0.0*cm;
  //G4double JawX2Pos_z = 54.2*cm;
  JawX2_phys = new G4PVPlacement(rotatejaw2X,
             G4ThreeVector(JawX2Pos_x,JawX2Pos_y,JawX2Pos_z),
             "JawX2",JawX2_log,experimentalHall_phys,false,0);


  //----------------Light Field Reticle

  G4double reticleDim_x = 15.0*cm;
  G4double reticleDim_y = 15.0*cm;
  G4double reticleDim_z = 0.00508*cm;
  G4Box* reticle = new G4Box("reticle",reticleDim_x,
                                  reticleDim_y,reticleDim_z);
  reticle_log = new G4LogicalVolume(reticle,
                                             Mylar,"reticle_log",0,0,0);

  G4double reticlePos_x = 0.0*cm;
  G4double reticlePos_y = 0.0*cm;
  G4double reticlePos_z = 50.5*cm;
  reticle_phys = new G4PVPlacement(0,
             G4ThreeVector(reticlePos_x,reticlePos_y,reticlePos_z),
             "reticle",reticle_log,experimentalHall_phys,false,0);

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
   Phantom_log->SetVisAttributes(simpleH2OVisAtt);

   G4VisAttributes* simpleTungstenWVisAtt= new G4VisAttributes(magenta);
   simpleTungstenWVisAtt->SetVisibility(true);
   simpleTungstenWVisAtt->SetForceWireframe(true);
   collim_log->SetVisAttributes(simpleTungstenWVisAtt);
   CylMinusCone_log->SetVisAttributes(simpleTungstenWVisAtt);
   UpperCollimator_log->SetVisAttributes(simpleTungstenWVisAtt);


   G4VisAttributes* simpleTungstenSVisAtt= new G4VisAttributes(magenta);
   simpleTungstenSVisAtt->SetVisibility(true);
   simpleTungstenSVisAtt->SetForceSolid(true);
   targetA_log->SetVisAttributes(simpleTungstenSVisAtt);
   JawY1_log->SetVisAttributes(simpleTungstenSVisAtt);
   JawY2_log->SetVisAttributes(simpleTungstenSVisAtt);
   JawX1_log->SetVisAttributes(simpleTungstenSVisAtt);
   JawX2_log->SetVisAttributes(simpleTungstenSVisAtt);

   G4VisAttributes* simpleCopperSVisAtt= new G4VisAttributes(cyan);
   simpleCopperSVisAtt->SetVisibility(true);
   simpleCopperSVisAtt->SetForceSolid(true);
   targetB_log->SetVisAttributes(simpleCopperSVisAtt);
   layer1_log->SetVisAttributes(simpleCopperSVisAtt);
   layer2_log->SetVisAttributes(simpleCopperSVisAtt);
   layer3_log->SetVisAttributes(simpleCopperSVisAtt);
   layer4_log->SetVisAttributes(simpleCopperSVisAtt);
   layer5_log->SetVisAttributes(simpleCopperSVisAtt);
   layer6_log->SetVisAttributes(simpleCopperSVisAtt);
   layer7_log->SetVisAttributes(simpleCopperSVisAtt);
   layer8_log->SetVisAttributes(simpleCopperSVisAtt);
   layer9_log->SetVisAttributes(simpleCopperSVisAtt);
   layer10_log->SetVisAttributes(simpleCopperSVisAtt);
   layer11_log->SetVisAttributes(simpleCopperSVisAtt);
   layer12_log->SetVisAttributes(simpleCopperSVisAtt);
   layer13_log->SetVisAttributes(simpleCopperSVisAtt);
   layer14_log->SetVisAttributes(simpleCopperSVisAtt);
   layer15_log->SetVisAttributes(simpleCopperSVisAtt);
   layer16_log->SetVisAttributes(simpleCopperSVisAtt);
   layer17_log->SetVisAttributes(simpleCopperSVisAtt);
   layer18_log->SetVisAttributes(simpleCopperSVisAtt);
   layer19_log->SetVisAttributes(simpleCopperSVisAtt);
   layer21_log->SetVisAttributes(simpleCopperSVisAtt);

   G4VisAttributes* simpleBerylliumSVisAtt= new G4VisAttributes(red);
   simpleBerylliumSVisAtt->SetVisibility(true);
   simpleBerylliumSVisAtt->SetForceSolid(true);
   windowUp_log->SetVisAttributes(simpleBerylliumSVisAtt);
   windowLow_log->SetVisAttributes(simpleBerylliumSVisAtt);

   G4VisAttributes* simpleMylarVisAtt= new G4VisAttributes(grey);
   simpleMylarVisAtt->SetVisibility(true);
   simpleMylarVisAtt->SetForceSolid(true);
   Mirror_log->SetVisAttributes(simpleBerylliumSVisAtt);
   reticle_log->SetVisAttributes(simpleMylarVisAtt);


   G4VisAttributes* simpleKaptonVisAtt= new G4VisAttributes(yellow);
   simpleKaptonVisAtt->SetVisibility(true);
   simpleKaptonVisAtt->SetForceSolid(true);
   Window_log->SetVisAttributes(simpleKaptonVisAtt);
   SignalPlate_log->SetVisAttributes(simpleKaptonVisAtt);


   G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(white);
   simpleWorldVisAtt->SetVisibility(false);
   experimentalHall_log ->SetVisAttributes(simpleWorldVisAtt);

   //se non voglio vedere phantom, jaws e mirror ...:

   //Phantom_log->SetVisAttributes(simpleWorldVisAtt);
   //JawY1_log->SetVisAttributes(simpleWorldVisAtt);
   //JawY2_log->SetVisAttributes(simpleWorldVisAtt);
   //JawX1_log->SetVisAttributes(simpleWorldVisAtt);
   //JawX2_log->SetVisAttributes(simpleWorldVisAtt);
   //Mirror_log->SetVisAttributes(simpleWorldVisAtt);
   //CylMinusCone_log->SetVisAttributes(simpleWorldVisAtt);
   //UpperCollimator_log->SetVisAttributes(simpleWorldVisAtt);
   //targetA_log->SetVisAttributes(simpleWorldVisAtt);
   //targetB_log->SetVisAttributes(simpleWorldVisAtt);
   //vacuumBlock_log->SetVisAttributes(simpleWorldVisAtt);
   //UpperCollimator_log->SetVisAttributes(simpleWorldVisAtt);
   //windowUp_log->SetVisAttributes(simpleWorldVisAtt);
   //windowLow_log->SetVisAttributes(simpleWorldVisAtt);
   //Window_log->SetVisAttributes(simpleWorldVisAtt);
   //SignalPlate_log->SetVisAttributes(simpleWorldVisAtt);
   //reticle_log->SetVisAttributes(simpleWorldVisAtt);


   //layer1_log->SetVisAttributes(simpleWorldVisAtt);
   //layer2_log->SetVisAttributes(simpleWorldVisAtt);
   //layer3_log->SetVisAttributes(simpleWorldVisAtt);
   //layer4_log->SetVisAttributes(simpleWorldVisAtt);
   //layer5_log->SetVisAttributes(simpleWorldVisAtt);
   //layer6_log->SetVisAttributes(simpleWorldVisAtt);
   //layer7_log->SetVisAttributes(simpleWorldVisAtt);
   //layer8_log->SetVisAttributes(simpleWorldVisAtt);
   //layer9_log->SetVisAttributes(simpleWorldVisAtt);
   //layer10_log->SetVisAttributes(simpleWorldVisAtt);
   //layer11_log->SetVisAttributes(simpleWorldVisAtt);
   //layer12_log->SetVisAttributes(simpleWorldVisAtt);
   //layer13_log->SetVisAttributes(simpleWorldVisAtt);
   //layer14_log->SetVisAttributes(simpleWorldVisAtt);
   //layer15_log->SetVisAttributes(simpleWorldVisAtt);
   //layer16_log->SetVisAttributes(simpleWorldVisAtt);
   //layer17_log->SetVisAttributes(simpleWorldVisAtt);
   //layer18_log->SetVisAttributes(simpleWorldVisAtt);
   //layer19_log->SetVisAttributes(simpleWorldVisAtt);
   //layer21_log->SetVisAttributes(simpleWorldVisAtt);
   //----------------------------------------------------------

   ConstructSensitiveDetector();
   return experimentalHall_phys;
}

void MedLinacDetectorConstruction::PrintParameters()
{
  G4cout <<"jaws1 x position "<< fieldX1/cm << " cm "<<G4endl ;
  G4cout <<"jaws2 x position "<< fieldX2/cm << " cm "<<G4endl ; 
  G4cout <<"jaws1 y position "<< fieldY1/cm << " cm "<<G4endl ;
  G4cout <<"jaws2 y position "<< fieldY2/cm << " cm "<<G4endl ; 
    
}



void MedLinacDetectorConstruction::SetJawX1Pos_x (G4double val)
{
  fieldX1 = val;
  G4cout <<"==============================DetectorConstruction JawX1 "<< fieldX1/cm<<"cm"<<G4endl;
}


void MedLinacDetectorConstruction::SetJawX2Pos_x (G4double val)
{
  fieldX2 = val;
  G4cout <<"==============================DetectorConstruction JawX2  "<< fieldX2/cm<<"cm"<<G4endl;
}

void MedLinacDetectorConstruction::SetJawY1Pos_y (G4double val)
{
  fieldY1 = val;
  G4cout <<"==============================DetectorConstruction JawY1 "<< fieldY1/cm<<"cm"<<G4endl;
}


void MedLinacDetectorConstruction::SetJawY2Pos_y (G4double val)
{
  fieldY2 = val;
  G4cout <<"==============================DetectorConstruction JawY2  "<< fieldY2/cm<<"cm"<<G4endl;
}


void MedLinacDetectorConstruction::UpdateGeometry()
{ 
 //delete experimentalHall_log;
 //delete experimentalHall_phys;

  G4cout <<"@@@@@@@@@@@@@@@@@@@@@UpdateGeometry "<< G4endl;

  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructGeom());
}



void  MedLinacDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  if(!phantomSD)
  {
    //G4double phantomDim_x = 15.*cm;
    //G4double phantomDim_y = 15.*cm;
    //G4double phantomDim_z = 15.*cm;
    phantomSD = new MedLinacPhantomSD(sensitiveDetectorName,numberOfVoxelsAlongX,numberOfVoxelsAlongY,numberOfVoxelsAlongZ);
    G4String ROGeometryName = "PhantomROGeometry";
    phantomROGeometry = new MedLinacPhantomROGeometry(ROGeometryName, phantomDim_x,phantomDim_y,phantomDim_z,
						  numberOfVoxelsAlongX,numberOfVoxelsAlongY,numberOfVoxelsAlongZ);

    //G4cout <<" numberOfVoxelsAlongX in det constr---------------------"<< numberOfVoxelsAlongX << G4endl;

    //G4cout <<"sensitiveDetectorName in det constr---------------------"<< sensitiveDetectorName << G4endl;

 G4cout <<"__________phantomDim_x---------------------"<< phantomDim_x/cm << G4endl;
 G4cout <<"__________numberOfVoxelsAlongX---------------------"<<numberOfVoxelsAlongX  << G4endl;
    phantomROGeometry->BuildROGeometry();
    phantomSD->SetROgeometry(phantomROGeometry);
    pSDManager->AddNewDetector(phantomSD);
    Phantom_log->SetSensitiveDetector(phantomSD);
  }
}














