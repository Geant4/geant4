//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: MedLinacHead.cc,v 1.3 2006/06/29 16:04:25 gunter Exp $
//
// Code developed by: M. Piergentili
//
//
//------------------------------ beam line along z axis---------------------

#include "MedLinacDetectorMessenger.hh"
#include "MedLinacDetectorConstruction.hh"
#include "MedLinacPhantomROGeometry.hh"
#include "MedLinacPhantomSD.hh"
#include "MedLinacVGeometryComponent.hh"
#include "MedLinacHead.hh"

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
#include "G4Region.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

MedLinacHead::MedLinacHead()
  : windowUp_log(0), UpperCollimator_log(0),
    collim_log(0), tracker_log(0),
    CylMinusCone_log(0),windowLow_log(0),
    Window_log(0),SignalPlate_log(0),
    Mirror_log(0),reticle_log(0),
    windowUp_phys(0), UpperCollimator_phys(0),
    CylMinusCone_phys(0), windowLow_phys(0),
    Window1_phys(0), SignalPlate1_phys(0),
    SignalPlate2_phys(0), Window2_phys(0),
    SignalPlate3_phys(0), SignalPlate4_phys(0),
    Window3_phys(0),Mirror_phys(0),
    reticle_phys(0), aLowerCollRegion(0),aUpperCollRegion(0)
   
{ 
  aLowerCollRegion = new G4Region("PrimaryCollimatorLow");
  aUpperCollRegion = new G4Region("PrimaryCollimatorUp");
}

MedLinacHead::~MedLinacHead()
{;}


void MedLinacHead::ConstructComponent(G4VPhysicalVolume* world,G4VPhysicalVolume* vacuumBlock)
{


  //G4GeometryManager::GetInstance()->OpenGeometry();
  //G4PhysicalVolumeStore::GetInstance()->Clean();
  //G4LogicalVolumeStore::GetInstance()->Clean();
  //G4SolidStore::GetInstance()->Clean();
  //G4RegionStore::GetInstance()->Clean(); // -> Segmentation fault (core dumped)


  //    materials

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4int ncomponents;
  //G4int natoms;
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

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air",density, ncomponents=2);
  Air->AddElement(elN, fractionmass=70*perCent);
  Air->AddElement(elO, fractionmass=30*perCent);

  //density = 2.702*g/cm3;
  //a = 26.981539*g/mole;
  //G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  //density = 8.960*g/cm3;
  //a = 63.55*g/mole;
  //G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

  //a = 47.88*g/mole;
  //density = 4.50*g/cm3;
  //Titanium = new G4Material("titanium" ,z = 22.,a,density);

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



  //------------------------------ beam line along z axis-------SSD=100cm-----------

 //---------rotation matrix first collimator and filter--------

  G4RotationMatrix*  rotateMatrix=new G4RotationMatrix();
  rotateMatrix->rotateX(180.0*deg);

  //---------rotation Mirror--------

  G4RotationMatrix*  rotateMirror=new G4RotationMatrix();
  rotateMirror->rotateY(35.0*deg);

  //---------------colors----------

  //G4Colour  white   (1.0, 1.0, 1.0);
  G4Colour  grey    (0.5, 0.5, 0.5);
  //G4Colour  lgrey   (.75, .75, .75);
  G4Colour  red     (1.0, 0.0, 0.0);
  //G4Colour  blue    (0.0, 0.0, 1.0);
  G4Colour  cyan    (0.0, 1.0, 1.0);
  G4Colour  magenta (1.0, 0.0, 1.0); 
  G4Colour  yellow  (1.0, 1.0, 0.0);
  //G4Colour  lblue   (0.20, .50, .85);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)
  //------------------------------ beam line along z axis---------------------  
  //--------------------window upper-------------------

  G4double windowUpDim_x = 4.*cm;
  G4double windowUpDim_y = 4.*cm;
  G4double windowUpDim_z = 0.0125*cm;
  G4Box* windowUp_box = new G4Box("windowUp_box",windowUpDim_x,windowUpDim_y,windowUpDim_z);
  windowUp_log = new G4LogicalVolume(windowUp_box,Be,"windowUp_log",0,0,0);
  G4double windowUpPos_x = 0.0*m;
  G4double windowUpPos_y = 0.0*m;
  G4double windowUpPos_z = 123.5225*cm;
  windowUp_phys = new G4PVPlacement(0,
            G4ThreeVector(windowUpPos_x,windowUpPos_y,windowUpPos_z),
            "windowUp",windowUp_log,world,false,0);


  //-------------------- the first collimator upper----------------


  G4double innerRadiusOfTheTubeEx = 1.0*cm;
  G4double outerRadiusOfTheTubeEx = 8.*cm;
  G4double hightOfTheTubeEx = 4.0*cm;
  G4double startAngleOfTheTubeEx = 0.*deg;
  G4double spanningAngleOfTheTubeEx = 360.*deg;
  G4Tubs* UpperCollimator = new G4Tubs("UpperCollimator",innerRadiusOfTheTubeEx,
                                    outerRadiusOfTheTubeEx,hightOfTheTubeEx,
				    startAngleOfTheTubeEx,spanningAngleOfTheTubeEx);
  UpperCollimator_log = new G4LogicalVolume(UpperCollimator,W,"UpperCollimator_log",0,0,0);

  G4double UpperCollimatorPosX = 0.*cm;
  G4double UpperCollimatorPosY = 0.*cm;
  G4double UpperCollimatorPosZ = 2.645*cm;

  UpperCollimator_phys = new G4PVPlacement(0,
					   G4ThreeVector(UpperCollimatorPosX,UpperCollimatorPosY,
							 UpperCollimatorPosZ),"UpperCollimator",
					   UpperCollimator_log,vacuumBlock,false,0);
  //-------------------- the first collimator lower----------------

  G4double  pRmin1 = 0.*cm;
  //G4double  pRmax1 = 0.3352897*cm;
  G4double  pRmax1 = 0.5*cm;
  G4double  pRmin2 = 0.*cm;
  G4double  pRmax2 = 1.7658592*cm;
  G4double  hightOfTheCone =3.2*cm;
  G4double  startAngleOfTheCone = 0.*deg;
  G4double  spanningAngleOfTheCone = 360.*deg;

  G4Cons* collim_cone = new G4Cons("collim_cone",pRmin1,pRmax1,pRmin2,
				   pRmax2,hightOfTheCone,startAngleOfTheCone,
				   spanningAngleOfTheCone);
  collim_log = new G4LogicalVolume(collim_cone,Vacuum,"collim_log",0,0,0);


  G4double innerRadiusOfTheTube = 0.*cm;
  G4double outerRadiusOfTheTube = 8.*cm;
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
  G4double CminusCPos_z = -4.455*cm;
  CylMinusCone_phys = new G4PVPlacement(rotateMatrix,
					G4ThreeVector(CminusCPos_x,CminusCPos_y,CminusCPos_z),
					"CylMinusCone",CylMinusCone_log,vacuumBlock,false,0);
  
  
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
  G4double windowLowPos_z = 105.9873*cm;

  windowLow_phys = new G4PVPlacement(0,
            G4ThreeVector(windowLowPos_x,windowLowPos_y,windowLowPos_z),
            "windowLow",windowLow_log,world,false,0);

     
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

  G4double Window1Pos_z = 100.165*cm;
  Window1_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,Window1Pos_z),
           "Window",Window_log,world,false,0);

 //------Ion chamber---Signal Plate 1---------------------

  G4double SignalPlate1Pos_z = 99.927*cm;
  SignalPlate1_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate1Pos_z),
           "SignalPlate",SignalPlate_log,world,false,0);

 //------Ion chamber---Signal Plate 2---------------------

  G4double SignalPlate2Pos_z = 99.688*cm;
  SignalPlate2_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate2Pos_z),
           "SignalPlate",SignalPlate_log,world,false,0);

 //------Ion chamber---window 2---------------------

  G4double Window2Pos_z = 99.45*cm;
  Window2_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,Window2Pos_z),
           "Window",Window_log,world,false,0);

 //------Ion chamber---Signal Plate 3---------------------

  G4double SignalPlate3Pos_z = 99.212*cm;
  SignalPlate3_phys = new G4PVPlacement(0,
					G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate3Pos_z),
					"SignalPlate",SignalPlate_log,world,false,0);

 //------Ion chamber---Signal Plate 4---------------------

  G4double SignalPlate4Pos_z = 98.973*cm;
  SignalPlate4_phys = new G4PVPlacement(0,
					G4ThreeVector(IonChamberPos_x,IonChamberPos_y,SignalPlate4Pos_z),
					"SignalPlate",SignalPlate_log,world,false,0);

 //------Ion chamber---window 3---------------------

  G4double Window3Pos_z = 98.735*cm;
  Window3_phys = new G4PVPlacement(0,
           G4ThreeVector(IonChamberPos_x,IonChamberPos_y,Window3Pos_z),
           "Window",Window_log,world,false,0);


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
  G4double MirrorPos_z = 93.0*cm;
  Mirror_phys = new G4PVPlacement(rotateMirror,
           G4ThreeVector(MirrorPos_x,MirrorPos_y,MirrorPos_z),
           "Mirror",Mirror_log,world,false,0);

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
  G4double reticlePos_z = 65.5*cm;
  reticle_phys = new G4PVPlacement(0,
             G4ThreeVector(reticlePos_x,reticlePos_y,reticlePos_z),
             "reticle",reticle_log,world,false,0);
  
  // Cuts by Regions 


  //G4String regName[] = {"PrimaryCollimatorUp","PrimaryCollimatorLow"};
  //G4Region* aUpperCollRegion = new G4Region(regName[0]);
  //G4Region* aLowerCollRegion = new G4Region(regName[1]);
  //UpperCollimator_log->SetRegion(aUpperCollRegion);
  //CylMinusCone_log->SetRegion(aLowerCollRegion);
  //aUpperCollRegion->AddRootLogicalVolume(UpperCollimator_log);
  //aLowerCollRegion->AddRootLogicalVolume(CylMinusCone_log);

  //G4String regName = "PrimaryCollimatorLow";
  //G4Region* aLowerCollRegion = new G4Region(regName);
  CylMinusCone_log->SetRegion(aLowerCollRegion);
  aLowerCollRegion->AddRootLogicalVolume(CylMinusCone_log);

  //G4String regName1 = "PrimaryCollimatorUp";
  //G4Region* aUpperCollRegion = new G4Region(regName1);
  UpperCollimator_log->SetRegion(aUpperCollRegion);
  aUpperCollRegion->AddRootLogicalVolume(UpperCollimator_log);

//--------- Visualization attributes -------------------------------

   G4VisAttributes* simpleTungstenWVisAtt= new G4VisAttributes(magenta);
   simpleTungstenWVisAtt->SetVisibility(true);
   simpleTungstenWVisAtt->SetForceWireframe(true);
   collim_log->SetVisAttributes(simpleTungstenWVisAtt);
   CylMinusCone_log->SetVisAttributes(simpleTungstenWVisAtt);
   UpperCollimator_log->SetVisAttributes(simpleTungstenWVisAtt);

   G4VisAttributes* simpleCopperSVisAtt= new G4VisAttributes(cyan);
   simpleCopperSVisAtt->SetVisibility(true);
   simpleCopperSVisAtt->SetForceSolid(true);
   
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

   //if you want not to see  mirror ...:

   //G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(white);
   //simpleWorldVisAtt->SetVisibility(false);
   //Mirror_log->SetVisAttributes(simpleWorldVisAtt);
   //CylMinusCone_log->SetVisAttributes(simpleWorldVisAtt);
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

}

void MedLinacHead::DestroyComponent()
{;}
