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
/// \file electromagnetic/TestEm10/src/Em10DetectorConstruction.cc
/// \brief Implementation of the Em10DetectorConstruction class
//
//
// $Id: Em10DetectorConstruction.cc 73033 2013-08-15 09:24:45Z gcosmo $
//
//

#include "Em10DetectorConstruction.hh"
#include "Em10DetectorMessenger.hh"
#include "Em10CalorimeterSD.hh"
#include "Em10Materials.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4ProductionCuts.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10DetectorConstruction::Em10DetectorConstruction()
  :G4VUserDetectorConstruction(),
  fWorldChanged(false), fAbsorberMaterial(0), fGapMat(0), fSetUp("simpleALICE"),
  fWorldMaterial(0), fSolidWorld(0), fLogicWorld(0), fPhysicsWorld(0),
//   fSolidRadSlice(0), fLogicRadSlice(0),  fPhysicRadSlice(0),
   fSolidRadiator(0),  fLogicRadiator(0),   fPhysicsRadiator(0),
   fRadiatorMat(0), fPipe(false), fPipeField(false),
   fSolidAbsorber(0),  fLogicAbsorber(0),   fPhysicsAbsorber(0),
   fMagField(0),       fCalorimeterSD(0),   fRegGasDet(0),
   fRadRegion(0), fMat(0)
{
  fDetectorMessenger = new Em10DetectorMessenger(this);
  fMat               = new Em10Materials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10DetectorConstruction::~Em10DetectorConstruction()
{
  delete fDetectorMessenger;
  delete fMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Em10DetectorConstruction::Construct()
{
  return ConstructDetectorXTR();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Em10DetectorConstruction::ConstructDetectorXTR()
{
 // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  if( fSetUp == "simpleALICE" )
  {
    return SimpleSetUpALICE();
  }
  else if( fSetUp == "alice06" )
  {
    return SetUpALICE06();
  }
  else if( fSetUp == "bari05" )
  {
    return SetUpBari05();
  }
  else if( fSetUp == "harris73" )
  {
    return SetUpHarris73();
  }
  else if( fSetUp == "watase86" )
  {
    return SetUpWatase86();
  }
  else if( fSetUp == "barr90" )
  {
    return SetUpBarr90();
  }
  else
  {
    G4cout <<
    "Experimental setup is unsupported. Check /XTRdetector/setup " <<G4endl;
    G4cout<<"Run default: barr90 "<<G4endl;
    return SetUpBarr90();

    //  return 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Simplified setup for ALICE XTR test beam (~2004). 
// Runs by : TestEm10 salice.mac

G4VPhysicalVolume* Em10DetectorConstruction::SimpleSetUpALICE()
{
  fWorldSizeZ = 400.*cm;
  fWorldSizeR = 20.*cm;

  // Radiator and detector parameters

  fRadThickness = 0.020*mm;
  fGasGap       = 0.250*mm;  
  foilGasRatio  = fRadThickness/(fRadThickness+fGasGap);

  fFoilNumber   = 220;

  fAbsorberThickness = 38.3*mm;

  fAbsorberRadius   = 100.*mm;
  fAbsorberZ        = 136.*cm;

  fWindowThick    = 51.0*micrometer;
  fElectrodeThick = 10.0*micrometer;
  fGapThick       =  10.0*cm;

  fDetThickness =  40.0*mm;
  fDetLength    = 200.0*cm;
  fDetGap       =   0.01*mm;

  fStartR       = 40*cm;
  fStartZ       = 100.0*mm;

  fModuleNumber = 1;

  // Preparation of mixed radiator material

  G4Material* Mylar = fMat->GetMaterial("Mylar");
  G4Material* Air   = fMat->GetMaterial("Air");
  G4Material* Al   = fMat->GetMaterial("Al");

  G4double foilDensity =  1.39*g/cm3;
  // Mylar // 0.91*g/cm3;  // CH2 0.534*g/cm3; //Li
  G4double gasDensity  =  1.2928*mg/cm3;
  // Air // 1.977*mg/cm3; // CO2 0.178*mg/cm3; //He
  G4double totDensity  = foilDensity*foilGasRatio + 
                                             gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;
 
  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,
                                                  2);
  radiatorMat->AddMaterial( Mylar, fractionFoil );
  radiatorMat->AddMaterial( Air, fractionGas  );

  // default materials of the detector and TR radiator

  fRadiatorMat =  radiatorMat;
  fFoilMat     = Mylar; // CH2; // Kapton; // Mylar ; // Li ; // CH2 ;
  fGasMat      = Air; // CO2; // He; //
  
  fWindowMat    = Mylar;
  fElectrodeMat = Al;

  fAbsorberMaterial = fMat->GetMaterial("Xe15CO2");
 
  fGapMat          = fAbsorberMaterial;

  fWorldMaterial    = Air; // CO2;

  fSolidWorld = new G4Box("World", fWorldSizeR,fWorldSizeR,fWorldSizeZ/2.);
 
  fLogicWorld = new G4LogicalVolume(fSolidWorld,  fWorldMaterial,  "World");

  fPhysicsWorld = new G4PVPlacement(0, G4ThreeVector(), "World",
                                 fLogicWorld, 0,  false, 0);

  // TR radiator envelope

  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap + fDetGap;

  fRadZ = fStartZ + 0.5*fRadThick;

  fSolidRadiator = new G4Box("Radiator",1.1*fAbsorberRadius ,
                              1.1*fAbsorberRadius,  0.5*fRadThick );

  fLogicRadiator = new G4LogicalVolume(fSolidRadiator, fRadiatorMat,
                                       "Radiator");
 
  fPhysicsRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,fRadZ),
                                     "Radiator", fLogicRadiator,
                                     fPhysicsWorld, false,        0 );

  // create region for window inside windowR for

  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);

  fWindowZ = fStartZ + fRadThick + fWindowThick/2. + 15.0*mm;

  //  G4Box* solidWindowR = new G4Box("WindowR",fAbsorberRadius+0.001,
  //                                        fAbsorberRadius+0.001,
  //                                        fWindowThick/2.+0.001  );

  //  G4LogicalVolume* logicWindowR = new G4LogicalVolume(solidWindowR,
  //                                   fWorldMaterial, "WindowR");
 
  //  G4VPhysicalVolume*    physiWindowR = new G4PVPlacement(0,
  //                      G4ThreeVector(0.,0.,fWindowZ),
  //                      "WindowR",logicWindowR,fPhysicsWorld,false,0);
  // window

  //  G4Box* solidWindow = new G4Box("Window",fAbsorberRadius,
  //                                 fAbsorberRadius, fWindowThick/2.);
 
  //  G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,
  //                                   fWindowMat, "Window");

  //  G4VPhysicalVolume* physiWindow = 
  //                         new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
  //                         "Window", logicWindow, physiWindowR, false, 0); 

  fGapZ = fWindowZ + fWindowThick/2. + fGapThick/2. + 0.01*mm;

  fElectrodeZ = fGapZ + fGapThick/2. + fElectrodeThick/2. + 0.01*mm;

  // Absorber

  fAbsorberZ = fElectrodeZ + fElectrodeThick/2. + 
                  fAbsorberThickness/2. + 0.01*mm;

  fSolidAbsorber = new G4Box("Absorber", fAbsorberRadius,
                                 fAbsorberRadius, fAbsorberThickness/2.);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber, fAbsorberMaterial,
                                                "Absorber");

  fPhysicsAbsorber = new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsorberZ),
                                       "Absorber", fLogicAbsorber,
                                        fPhysicsWorld,  false,  0);

  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

  // Sensitive Detectors: Absorber

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  return fPhysicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Setup for ALICE XTR test beam (~2004). With He beam-pipe
// Runs by : TestEm10 alice06.mac

G4VPhysicalVolume* Em10DetectorConstruction::SetUpALICE06()
{
  fWorldSizeZ = 600.*cm;
  fWorldSizeR = 22.*cm;

  // Radiator and detector parameters

  //fRadThickness = 0.01*mm;    // Gamma XTR (malz: 0.01)
  //fGasGap       = 0.19*mm;    // Gamma XTR (malz: 0.09)
  //fFoilNumber   = 240;        // Gamma XTR (malz: 480)

  fRadThickness = 0.020*mm;  // Reg1
  fGasGap       = 0.500*mm;  // Reg1
  fFoilNumber   = 120;       // Reg1

  //fRadThickness = 0.013*mm;  // Anton
  //fGasGap       = 0.060*mm;  // Anton
  //fFoilNumber   = 550;       // Anton

  // fRadThickness = 0.020*mm; // Reg2
  // fGasGap       = 0.250*mm; // Reg2 
  // fFoilNumber   = 220;      // Reg2

  foilGasRatio  = fRadThickness/(fRadThickness+fGasGap);

  fAbsorberThickness = 37.*mm; // 38.3*mm;

  fAbsorberRadius   = 100.*mm;
  fAbsorberZ        = 136.*cm;

  fPipeLength     = 160.0*cm;
  fMylarThick     = 20.0*micrometer;

  fWindowThick    = 51.0*micrometer;
  fElectrodeThick = 100.0*micrometer;
  fGapThick       =  10.0*cm;

  fDetThickness =  40.0*mm;
  fDetLength    = 200.0*cm;
  fDetGap       =   0.01*mm;

  fStartR       = 40*cm;
  fStartZ       = 100.0*mm;

  fModuleNumber = 1;

  // Preparation of mixed radiator material

  G4Material* Mylar = fMat->GetMaterial("Mylar");
  G4Material* Air   = fMat->GetMaterial("Air");
  G4Material* Al   = fMat->GetMaterial("Al");
  G4Material* CH2   = fMat->GetMaterial("CH2");
  G4Material* He   = fMat->GetMaterial("He");

  G4double foilDensity = CH2->GetDensity();
  G4double gasDensity  = Air->GetDensity();  
  G4double totDensity  = foilDensity*foilGasRatio + 
                                              gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  1.0 - fractionFoil;
// gasDensity*(1.0-foilGasRatio)/totDensity ;  
    
  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,
                                                  2);
  radiatorMat->AddMaterial( CH2, fractionFoil );
  radiatorMat->AddMaterial( Air, fractionGas  );

  // default materials of the detector and TR radiator

  fRadiatorMat = radiatorMat;
  fFoilMat     = CH2; // Kapton; // Mylar ; // Li ; // CH2 ;
  fGasMat      = Air; // CO2; // He; //

  fWindowMat    = Mylar;
  fElectrodeMat = Al;

  fAbsorberMaterial = fMat->GetMaterial("Xe15CO2");

  // pipe material is assumed to be He + small admixture of air
  /* 
  foilGasRatio = 0.000001;
  foilDensity  = Air->GetDensity();
  gasDensity   = He->GetDensity();
  totDensity   = foilDensity*foilGasRatio + gasDensity*( 1.0 - foilGasRatio );

  fractionFoil =  foilDensity*foilGasRatio/totDensity;
  fractionGas  =  1.0 - fractionFoil;
  // gasDensity*(1.0 - foilGasRatio)/totDensity;  

  fPipeMat = new G4Material("pipeMat"  , totDensity,  2);
  fPipeMat->AddMaterial( Air, fractionFoil );
  fPipeMat->AddMaterial( He,  fractionGas  );
  */
  fPipeMat = He;

  fGapMat           = fAbsorberMaterial;

  fWorldMaterial    = Air;

  fSolidWorld = new G4Box("World", fWorldSizeR, fWorldSizeR, fWorldSizeZ/2.);

  fLogicWorld = new G4LogicalVolume(fSolidWorld,  fWorldMaterial,  "World");

  fPhysicsWorld = new G4PVPlacement(0, G4ThreeVector(), "World",
                                 fLogicWorld, 0,  false, 0);

  // TR radiator envelope

  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap + fDetGap;

  fRadZ = fStartZ + 0.5*fRadThick;

  // fRadZ = -fRadThick/2. - fElectrodeThick;
  // if ( fabs(pipe) > 1.e-15 ) fRadZ -= ( fPipeLength/2. + pipeDist );

  fSolidRadiator = new G4Box("Radiator",1.1*fAbsorberRadius ,
                              1.1*fAbsorberRadius,  0.5*fRadThick );

  fLogicRadiator = new G4LogicalVolume(fSolidRadiator, fRadiatorMat,
                                       "Radiator");

  fPhysicsRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,fRadZ),
                                     "Radiator", fLogicRadiator,
                                     fPhysicsWorld, false,        0 );

  // create region for radiator

  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);

  // Drift Electrode on both sides of Radiator:

  G4double zElectrode1 = fRadZ - fRadThick/2. - fElectrodeThick/2.;
  G4double zElectrode2 = fRadZ + fRadThick/2. + fElectrodeThick/2.;
  /*
  G4Box* solidElectrode = new G4Box("Electrode",fAbsorberRadius*1.1,
                                                fAbsorberRadius*1.1,
                                                fElectrodeThick/2.);

  G4LogicalVolume* logicElectrode = new G4LogicalVolume(solidElectrode,
                                                        fElectrodeMat,
                                                        "Electrode");

  G4VPhysicalVolume*    physiElectrode1 = new G4PVPlacement(0,
                                       G4ThreeVector(0.,0.,zElectrode1),
                                      "Electrode1",logicElectrode,
                                       fPhysicsWorld,false,0);

  G4VPhysicalVolume*    physiElectrode2 = new G4PVPlacement(0,
                                       G4ThreeVector(0.,0.,zElectrode2),
                                      "Electrode1",logicElectrode,
                                       fPhysicsWorld,false,0);
  */
  G4cout<<"zElectrode1 = "<<zElectrode1/mm<<" mm"<<G4endl;
  G4cout<<"zElectrode2 = "<<zElectrode2/mm<<" mm"<<G4endl;
  G4cout<<"fElectrodeThick = "<<fElectrodeThick/mm<<" mm"<<G4endl<<G4endl;

  // Helium Pipe:

  //Distance between pipe and radiator / absorber
  G4double pipeDist      = 1.*cm;
  G4double fieldStrength = 1.0*tesla;  // 0.01*tesla; // field strength in pipe
  G4double alphaB        = 90.*degree;
  fPipe     =  true;   // 0.;  //  use helium pipe is setup

  fPipeField     =  true;   // field in helium pipe used?

  G4double zPipe = zElectrode2 + fElectrodeThick/2. +
                   pipeDist/2. + fPipeLength/2.;

  if ( fPipe )
  {

    G4Box* solidPipe = new G4Box("Pipe",fAbsorberRadius*0.5,
                                 fAbsorberRadius*0.5,
                                 fPipeLength/2. );

    G4LogicalVolume* logicPipe = new G4LogicalVolume(solidPipe,
                                                     fPipeMat, //fWorldMaterial
                                                     "Pipe");

    //    G4VPhysicalVolume*    physiPipe = new G4PVPlacement(0,
    //                                 G4ThreeVector(0., 0., zPipe),
    //                                "Pipe1",logicPipe,
    //                                  fPhysicsWorld,false,0);

    G4cout<<"zPipe = "<<zPipe/mm<<" mm"<<G4endl;
    G4cout<<"fPipeLength = "<<fPipeLength/mm<<" mm"<<G4endl<<G4endl;

    // magnetic field in Pipe:

    if ( fPipeField )
    {
      if( fMagField ) delete fMagField; //delete the existing mag field

       fMagField =
           new G4UniformMagField(G4ThreeVector(fieldStrength*std::sin(alphaB), 
                                 0., fieldStrength*std::cos(alphaB)));
      // fMagField = new G4UniformMagField(G4ThreeVector(fieldStrength,0.,0.));
      // fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldStrength));
      G4FieldManager* fieldMgr = new G4FieldManager(fMagField);
      fieldMgr->SetDetectorField(fMagField);
      fieldMgr->CreateChordFinder(fMagField);
      logicPipe->SetFieldManager(fieldMgr, true);
    }

  }
  else   G4cout<<"No Helium pipe is used"<<G4endl<<G4endl;

  // Mylar Foil on both sides of helium pipe:

  G4double zMylar1 = zPipe - fPipeLength/2. - fMylarThick/2. - 0.001*mm;
  G4double zMylar2 = zPipe + fPipeLength/2. + fMylarThick/2. + 0.001*mm;

  //  G4Box* solidMylar = new G4Box("MylarB",fAbsorberRadius*0.6,
  //                              fAbsorberRadius*0.6,
  //                               fMylarThick/2.);

  //  G4LogicalVolume* logicMylar = new G4LogicalVolume(solidMylar,
  //                                                  fWindowMat,
  //                                                  "MylarL");

  if ( fPipe )
  {

    //    G4VPhysicalVolume* physiMylar1 = new G4PVPlacement(0,
    //                         G4ThreeVector( 0., 0., zMylar1),
    //                          "Mylar1", logicMylar, fPhysicsWorld,
    //                                      false, 0);

    //  G4VPhysicalVolume* physiMylar2 = new G4PVPlacement(0,
    //                             G4ThreeVector(0., 0., zMylar2),
    //                            "Mylar2", logicMylar, fPhysicsWorld,
    //                               false, 0);

      G4cout<<"zMylar1 = "<<zMylar1/mm<<" mm"<<G4endl;
      G4cout<<"zMylar2 = "<<zMylar2/mm<<" mm"<<G4endl;
      G4cout<<"fMylarThick = "<<fMylarThick/mm<<" mm"<<G4endl<<G4endl;
  }

  // Mylar Foil on Chamber:

  G4double zMylar = zElectrode2 + fElectrodeThick/2. + fMylarThick/2. + 1.0*mm;

  // if ( fPipe )
  {
    zMylar += ( fPipeLength + pipeDist );
  }
  //  G4VPhysicalVolume*    physiMylar = new G4PVPlacement(0,
  //                       G4ThreeVector(0., 0., zMylar),
  //                      "Mylar",logicMylar,fPhysicsWorld,false,0);

  G4cout<<"zMylar = "<<zMylar/mm<<" mm"<<G4endl;
  G4cout<<"fMylarThick = "<<fMylarThick/mm<<" mm"<<G4endl<<G4endl;

  // Absorber

  fAbsorberZ = zMylar + fMylarThick + fAbsorberThickness/2.;

  fSolidAbsorber = new G4Box("Absorber",
                             fAbsorberRadius,
                             // fAbsorberRadius,
                             // 10.*mm,
                             10.*mm,
                             fAbsorberThickness/2.);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber, fAbsorberMaterial,
                                                "Absorber");

  fPhysicsAbsorber = new G4PVPlacement(0,
                         G4ThreeVector(0., 0., fAbsorberZ),
                                       "Absorber", fLogicAbsorber,
                                        fPhysicsWorld,  false,  0);

  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

  // Sensitive Detectors: Absorber

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  return fPhysicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Setup for Bari INFN XTR test beam (~2004) at CERN. With He beam-pipe
// M. Brigida et al, NIM A550 (2005) 157-168
// Runs by : TestEm10 bari05.mac

G4VPhysicalVolume* Em10DetectorConstruction::SetUpBari05()
{
  fWorldSizeZ = 600.*cm;
  fWorldSizeR = 22.*cm;

  // Radiator and detector parameters

  //fRadThickness = 0.01*mm;    // Gamma XTR (malz: 0.01)
  //fGasGap       = 0.19*mm;    // Gamma XTR (malz: 0.09)
  //fFoilNumber   = 240;        // Gamma XTR (malz: 480)

  //fRadThickness = 0.020*mm;  // Reg1
  //fGasGap       = 0.500*mm;  // Reg1
  //fFoilNumber   = 120;       // Reg1

  //fRadThickness = 0.013*mm;  // Anton
  //fGasGap       = 0.230*mm;  // Anton
  //fFoilNumber   = 550;       // Anton

  fRadThickness = 0.0055*mm; // Reg2
  fGasGap       = 0.23*mm; // Reg2
  fFoilNumber   = 191;      // Reg2

  foilGasRatio  = fRadThickness/(fRadThickness+fGasGap);

  fAbsorberThickness = 0.4*mm;

  fAbsorberRadius   = 100.*mm;
  fAbsorberZ        = 136.*cm;

  fPipeLength = 50.0*cm;
  fMylarThick     = 20.0*micrometer;

  fWindowThick    = 51.0*micrometer;
  fElectrodeThick = 100.0*micrometer;
  fGapThick       =  10.0*cm;

  fDetThickness =  40.0*mm;
  fDetLength    = 200.0*cm;
  fDetGap       =   0.01*mm;

  fStartR       = 40*cm;
  fStartZ       = 100.0*mm;

  fModuleNumber = 1;

  // Preparation of mixed radiator material

  G4Material* Mylar = fMat->GetMaterial("Mylar");
  G4Material* Air   = fMat->GetMaterial("Air");
  G4Material* Al   = fMat->GetMaterial("Al");
  G4Material* CH2   = fMat->GetMaterial("CH2");
  G4Material* He   = fMat->GetMaterial("He");

  G4double foilDensity =  0.91*g/cm3;  
  // CH2 1.39*g/cm3; // Mylar //  0.534*g/cm3; //Li
  G4double gasDensity  =  1.2928*mg/cm3;
  // Air // 1.977*mg/cm3; // CO2 0.178*mg/cm3; // He
  G4double totDensity  = foilDensity*foilGasRatio + 
                                           gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;

  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,
                                                  2);
  radiatorMat->AddMaterial( CH2, fractionFoil );
  radiatorMat->AddMaterial( Air, fractionGas  );

  // default materials of the detector and TR radiator

  fRadiatorMat = radiatorMat;
  fFoilMat     = CH2; // Kapton; // Mylar ; // Li ; // CH2 ;
  fGasMat      = Air; // CO2; // He; //

  fWindowMat    = Mylar;
  fElectrodeMat = Al;

  fAbsorberMaterial = fMat->GetMaterial("Si");

  // pipe material is assumed to be He + small admixture of air

  foilGasRatio = 0.99999;
  foilDensity  = 1.2928*mg/cm3; // Air
  gasDensity   = 0.178*mg/cm3; // He
  totDensity   = foilDensity*foilGasRatio + gasDensity*(1.0-foilGasRatio);

  fractionFoil =  foilDensity*foilGasRatio/totDensity;
  fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;

  fPipeMat = new G4Material("pipeMat"  , totDensity,  2);
  fPipeMat->AddMaterial( Air, fractionFoil );
  fPipeMat->AddMaterial( He,  fractionGas  );

  // fPipeMat = He;

  fGapMat           = fAbsorberMaterial;

  fWorldMaterial    = Air;

  fSolidWorld = new G4Box("World", fWorldSizeR,fWorldSizeR,fWorldSizeZ/2.);

  fLogicWorld = new G4LogicalVolume(fSolidWorld,  fWorldMaterial,  "World");

  fPhysicsWorld = new G4PVPlacement(0, G4ThreeVector(), "World",
                                 fLogicWorld, 0,  false, 0);

  // TR radiator envelope

  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap + fDetGap;

  fRadZ = fStartZ + 0.5*fRadThick;
  // fRadZ = -fRadThick/2. - fElectrodeThick;
  // if ( fabs(pipe) > 1.e-15 ) fRadZ -= ( fPipeLength/2. + pipeDist );

  fSolidRadiator = new G4Box("Radiator",1.1*fAbsorberRadius ,
                              1.1*fAbsorberRadius,  0.5*fRadThick );

  fLogicRadiator = new G4LogicalVolume(fSolidRadiator, fRadiatorMat,
                                       "Radiator");

  fPhysicsRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,fRadZ),
                                     "Radiator", fLogicRadiator,
                                     fPhysicsWorld, false,        0 );

  // create region for radiator

  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);

  // Drift Electrode on both sides of Radiator:

  //  G4Box* solidElectrode = new G4Box("Electrode",fAbsorberRadius*1.1,
  //                                            fAbsorberRadius*1.1,
  //                                             fElectrodeThick/2.);

  //  G4LogicalVolume* logicElectrode = new G4LogicalVolume(solidElectrode,
  //                                                       fElectrodeMat,
  //                                                        "Electrode");

  G4double zElectrode1 = fRadZ - fRadThick/2. - fElectrodeThick/2.;
  G4double zElectrode2 = fRadZ + fRadThick/2. + fElectrodeThick/2.;

  //  G4VPhysicalVolume*    physiElectrode1 = new G4PVPlacement(0,
  //                                       G4ThreeVector(0.,0.,zElectrode1),
  //                                     "Electrode1",logicElectrode,
  //                                      fPhysicsWorld,false,0);

  // G4VPhysicalVolume*    physiElectrode2 = new G4PVPlacement(0,
  //                                      G4ThreeVector(0.,0.,zElectrode2),
  //                                    "Electrode1",logicElectrode,
  //                                     fPhysicsWorld,false,0);

  G4cout<<"zElectrode1 = "<<zElectrode1/mm<<" mm"<<G4endl;
  G4cout<<"zElectrode2 = "<<zElectrode2/mm<<" mm"<<G4endl;
  G4cout<<"fElectrodeThick = "<<fElectrodeThick/mm<<" mm"<<G4endl<<G4endl;

  // Helium Pipe:

  G4double pipe     = 1.0;   // use helium pipe is setup

  G4double pipeDist = 1.*cm;  //Distance between pipe and radiator / absorber

  G4double zPipe = zElectrode2 + fElectrodeThick/2. + 
                                                fPipeLength/2. + pipeDist/2.;

  // G4double field         = 1.0;   // field in helium pipe used?
  // G4double fieldStrength = 1.0*tesla;  // field strength in pipe

  if ( std::fabs(pipe) > 1.e-15 )
  {

    //    G4Box* solidPipe = new G4Box("Pipe",fAbsorberRadius*0.5,
    //                              fAbsorberRadius*0.5,
    //                              fPipeLength/2. );

    //    G4LogicalVolume* logicPipe = new G4LogicalVolume(solidPipe,
    //                                                  fPipeMat,
    //                                                  "Pipe");

    // magnetic field in Pipe:
    // if( fMagField ) delete fMagField; //delete the existing mag field
    // fMagField = new G4UniformMagField(G4ThreeVector(fieldStrength,0.,0.));
    // G4FieldManager* fieldMgr= new G4FieldManager(fMagField);
    // fieldMgr->SetDetectorField(fMagField);
    // fieldMgr->CreateChordFinder(fMagField);
    // if ( fabs(field) > 1.e-15 ) logicPipe->SetFieldManager(fieldMgr, true);

    //    G4VPhysicalVolume*    physiPipe = new G4PVPlacement(0,
    //                                   G4ThreeVector(0.,0.,zPipe),
    //                                  "Pipe1",logicPipe,
    //                                   fPhysicsWorld,false,0);

    G4cout<<"zPipe = "<<zPipe/mm<<" mm"<<G4endl;
    G4cout<<"fPipeLength = "<<fPipeLength/mm<<" mm"<<G4endl<<G4endl;

  }
  else   G4cout<<"No Helium pipe is used"<<G4endl<<G4endl;

  // Mylar Foil on both sides of helium pipe:

  G4double zMylar1 = zPipe - fPipeLength/2. - fMylarThick/2 - 0.01*mm;
  G4double zMylar2 = zPipe + fPipeLength/2. + fMylarThick/2 + 0.01*mm;

  //  G4Box* solidMylar = new G4Box("Mylar",fAbsorberRadius*0.6,
  //                              fAbsorberRadius*0.6,
  //                              fMylarThick/2.);

  //  G4LogicalVolume* logicMylar = new G4LogicalVolume(solidMylar,
  //                                                  fWindowMat,
  //                                                  "Mylar");

  if ( std::fabs(pipe) > 1.e-15 )
  {

    //    G4VPhysicalVolume* physiMylar1 = new G4PVPlacement(0,
    //                           G4ThreeVector( 0., 0., zMylar1),
    //                            "Mylar1", logicMylar, fPhysicsWorld,
    //                                        false, 0);

    //  G4VPhysicalVolume* physiMylar2 = new G4PVPlacement(0,
    //                             G4ThreeVector(0.,0.,zMylar2),
    //                             "Mylar2", logicMylar, fPhysicsWorld,
    //                                false, 0);

      G4cout<<"zMylar1 = "<<zMylar1/mm<<" mm"<<G4endl;
      G4cout<<"zMylar2 = "<<zMylar2/mm<<" mm"<<G4endl;
      G4cout<<"fMylarThick = "<<fMylarThick/mm<<" mm"<<G4endl<<G4endl;

  }

  // Mylar Foil on Chamber:

  G4double zMylar = zElectrode2 + fElectrodeThick/2. + fMylarThick/2. + 1.0*mm;

  if ( std::fabs(pipe) > 1.e-15 ) zMylar += ( fPipeLength + pipeDist );

  //  G4VPhysicalVolume*    physiMylar = new G4PVPlacement(0,
  //                       G4ThreeVector(0.,0.,zMylar),
  //                      "Mylar",logicMylar,fPhysicsWorld,false,0);

  G4cout<<"zMylar = "<<zMylar/mm<<" mm"<<G4endl;
  G4cout<<"fMylarThick = "<<fMylarThick/mm<<" mm"<<G4endl<<G4endl;

  // Absorber

  fAbsorberZ = zMylar + fMylarThick/2. + fAbsorberThickness/2.;

  fSolidAbsorber = new G4Box("Absorber",
                             // fAbsorberRadius, fAbsorberRadius,
                             10.*mm,10.*mm,
                              fAbsorberThickness/2.);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber, fAbsorberMaterial,
                                                "Absorber");

  fPhysicsAbsorber = new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsorberZ),
                                       "Absorber", fLogicAbsorber,
                                        fPhysicsWorld,  false,  0);

  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");  
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

  // Sensitive Detectors: Absorber

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  return fPhysicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Setuo from F. Harris et al NIM 107 (1973) 413-422 (fig.5b)

G4VPhysicalVolume* Em10DetectorConstruction::SetUpHarris73()
{
  fWorldSizeZ = 400.*cm;
  fWorldSizeR = 20.*cm;

  // Radiator and detector parameters

  fRadThickness = 0.0127*mm;
  fGasGap       = 0.762*mm;
  foilGasRatio  = fRadThickness/(fRadThickness+fGasGap);

  fFoilNumber   = 100;

  fAbsorberThickness = 15.0*mm;

  fAbsorberRadius   = 100.*mm;
  fAbsorberZ        = 136.*cm;

  fWindowThick    = 51.0*micrometer;
  fElectrodeThick = 10.0*micrometer;
  fGapThick       =  10.0*cm;

  fDetThickness =  40.0*mm;
  fDetLength    = 200.0*cm;
  fDetGap       =   0.01*mm;

  fStartR       = 40*cm;
  fStartZ       = 100.0*mm;

  fModuleNumber = 1;

  // Preparation of mixed radiator material

  G4Material* Mylar = fMat->GetMaterial("Mylar");
  G4Material* Air   = fMat->GetMaterial("Air");
  G4Material* Al   = fMat->GetMaterial("Al");

  G4double foilDensity =  1.39*g/cm3;
  // Mylar // 0.91*g/cm3;  // CH2 0.534*g/cm3; //Li
  G4double gasDensity  =  1.2928*mg/cm3;
  // Air // 1.977*mg/cm3; // CO2 0.178*mg/cm3; // He
 
  G4double totDensity  = foilDensity*foilGasRatio +
                                            gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;

  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,
                                                  2);
  radiatorMat->AddMaterial( Mylar, fractionFoil );
  radiatorMat->AddMaterial( Air, fractionGas  );

  // default materials of the detector and TR radiator

  fRadiatorMat =  radiatorMat;
  fFoilMat     = Mylar;
  fGasMat      = Air;

  fWindowMat    = Mylar;
  fElectrodeMat = Al;

  fAbsorberMaterial = fMat->GetMaterial("Kr7CH4");
 
  fGapMat          = fAbsorberMaterial;

  fWorldMaterial    = Air; // CO2;

  fSolidWorld = new G4Box("World", fWorldSizeR,fWorldSizeR,fWorldSizeZ/2.);
 
  fLogicWorld = new G4LogicalVolume(fSolidWorld,  fWorldMaterial,  "World");

  fPhysicsWorld = new G4PVPlacement(0, G4ThreeVector(), "World",
                                 fLogicWorld, 0,  false, 0);

  // TR radiator envelope

  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap + fDetGap;

  fRadZ = fStartZ + 0.5*fRadThick;

  fSolidRadiator = new G4Box("Radiator",1.1*fAbsorberRadius ,
                              1.1*fAbsorberRadius,  0.5*fRadThick );

  fLogicRadiator = new G4LogicalVolume(fSolidRadiator, fRadiatorMat,
                                       "Radiator");

  fPhysicsRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,fRadZ),
                                     "Radiator", fLogicRadiator,
                                     fPhysicsWorld, false,        0 );

  // create region for window inside windowR for

  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);
 
  fWindowZ = fStartZ + fRadThick + fWindowThick/2. + 15.0*mm;

  // G4Box* solidWindowR = new G4Box("WindowR",fAbsorberRadius+0.001,
  //                                        fAbsorberRadius+0.001,
  //                                        fWindowThick/2.+0.001  ); 

  //  G4LogicalVolume* logicWindowR = new G4LogicalVolume(solidWindowR,
  //                                    fWorldMaterial, "WindowR");
 
  //  G4VPhysicalVolume*    physiWindowR = new G4PVPlacement(0,
  //                       G4ThreeVector(0.,0.,fWindowZ),
  //                             "WindowR",logicWindowR,fPhysicsWorld,false,0);
  // window
 
  // G4Box* solidWindow = new G4Box("Window",fAbsorberRadius,
  //                                  fAbsorberRadius, fWindowThick/2.);

  //  G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,
  //                                   fWindowMat, "Window");

  //  G4VPhysicalVolume*    physiWindow = 
  //                        new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
  //                        "Window", logicWindow, physiWindowR, false, 0);

  fGapZ = fWindowZ + fWindowThick/2. + fGapThick/2. + 0.01*mm;

  fElectrodeZ = fGapZ + fGapThick/2. + fElectrodeThick/2. + 0.01*mm;

  // Absorber

  fAbsorberZ = fElectrodeZ + fElectrodeThick/2. +
                                             fAbsorberThickness/2. + 0.01*mm;

  fSolidAbsorber = new G4Box("Absorber", fAbsorberRadius,
                                 fAbsorberRadius, fAbsorberThickness/2.);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber, fAbsorberMaterial,
                                                "Absorber");

  fPhysicsAbsorber = new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsorberZ),
                                       "Absorber", fLogicAbsorber,
                                        fPhysicsWorld,  false,  0);

  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");  
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

  // Sensitive Detectors: Absorber

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  return fPhysicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Setuo from Y. Watase et al, NIM A248  (1986) 379-388 (fig.7; Li, e-, 2 Gev/c)

G4VPhysicalVolume* Em10DetectorConstruction::SetUpWatase86()
{
  fWorldSizeZ = 400.*cm;
  fWorldSizeR = 20.*cm;

  // Radiator and detector parameters

  fRadThickness = 0.04*mm;
  fGasGap       = 0.126*mm;
  foilGasRatio  = fRadThickness/(fRadThickness+fGasGap);

  fFoilNumber   = 300;

  fAbsorberThickness = 30.0*mm;

  fAbsorberRadius   = 100.*mm;
  fAbsorberZ        = 136.*cm;

  fWindowThick    = 51.0*micrometer;
  fElectrodeThick = 10.0*micrometer;
  fGapThick       =  10.0*cm;

  fDetThickness =  30.0*mm;
  fDetLength    = 200.0*cm;
  fDetGap       =   0.01*mm;

  fStartR       = 40*cm;
  fStartZ       = 100.0*mm;

  fModuleNumber = 1;

  // Preparation of mixed radiator material

  G4Material* Li = fMat->GetMaterial("Li");
  //  G4Material* Air   = fMat->GetMaterial("Air");
  G4Material* He   = fMat->GetMaterial("He");
  G4Material* Al   = fMat->GetMaterial("Al");
  G4Material* Mylar = fMat->GetMaterial("Mylar");

  G4double foilDensity = 0.534*g/cm3;
  //Li  1.39*g/cm3; // Mylar 0.91*g/cm3;  // CH2
  G4double gasDensity  = 0.178*mg/cm3;
  // He 1.2928*mg/cm3; // Air // 1.977*mg/cm3; // CO2
 
  G4double totDensity  = foilDensity*foilGasRatio + 
                                            gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;

  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,
                                                  2);
  radiatorMat->AddMaterial( Li, fractionFoil );
  radiatorMat->AddMaterial( He, fractionGas  );

  // default materials of the detector and TR radiator

  fRadiatorMat =  radiatorMat;
  fFoilMat     = Li;
  fGasMat      = He;  

  fWindowMat    = Mylar;
  fElectrodeMat = Al;

  fAbsorberMaterial = fMat->GetMaterial("Xe10CH4");
 
  fGapMat          = fAbsorberMaterial;

  fWorldMaterial    = He; // Air; // CO2 ;

  fSolidWorld = new G4Box("World", fWorldSizeR,fWorldSizeR,fWorldSizeZ/2.);

  fLogicWorld = new G4LogicalVolume(fSolidWorld,  fWorldMaterial,  "World");

  fPhysicsWorld = new G4PVPlacement(0, G4ThreeVector(), "World",
                                 fLogicWorld, 0,  false, 0);

  // TR radiator envelope

  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap + fDetGap;

  fRadZ = fStartZ + 0.5*fRadThick;

  fSolidRadiator = new G4Box("Radiator",1.1*fAbsorberRadius ,
                              1.1*fAbsorberRadius,  0.5*fRadThick );

  fLogicRadiator = new G4LogicalVolume(fSolidRadiator, fRadiatorMat,
                                       "Radiator");

  fPhysicsRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,fRadZ),
                                     "Radiator", fLogicRadiator,
                                     fPhysicsWorld, false,        0 );

  // create region for window inside windowR for

  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);
 
  fWindowZ = fStartZ + fRadThick + fWindowThick/2. + 15.0*mm;

  // G4Box* solidWindowR = new G4Box("WindowR",fAbsorberRadius+0.001,
  //                                         fAbsorberRadius+0.001,
  //                                         fWindowThick/2.+0.001  );

  // G4LogicalVolume* logicWindowR = new G4LogicalVolume(solidWindowR,
  //                                    fWorldMaterial, "WindowR");
 
  //  G4VPhysicalVolume*    physiWindowR = new G4PVPlacement(0,
  //                    G4ThreeVector(0.,0.,fWindowZ),
  //                          "WindowR",logicWindowR,fPhysicsWorld,false,0);
  // window
 
  // G4Box* solidWindow = new G4Box("Window",fAbsorberRadius,
  //                                 fAbsorberRadius, fWindowThick/2.);

  //  G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,
  //                                    fWindowMat, "Window");

  //  G4VPhysicalVolume*    physiWindow =
  //                        new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
  //                        "Window", logicWindow, physiWindowR, false, 0);

  fGapZ = fWindowZ + fWindowThick/2. + fGapThick/2. + 0.01*mm;

  fElectrodeZ = fGapZ + fGapThick/2. + fElectrodeThick/2. + 0.01*mm;

  // Absorber

  fAbsorberZ = fElectrodeZ + fElectrodeThick/2. +
                                              fAbsorberThickness/2. + 0.01*mm;

  fSolidAbsorber = new G4Box("Absorber", fAbsorberRadius,
                                 fAbsorberRadius, fAbsorberThickness/2.);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber, fAbsorberMaterial,
                                                "Absorber");

  fPhysicsAbsorber = new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsorberZ),
                                       "Absorber", fLogicAbsorber,
                                        fPhysicsWorld,  false,  0);

  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");  
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

  // Sensitive Detectors: Absorber

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  return fPhysicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Setuo from G.D. Barr et al NIM A294 (1990) 465-472 (fig.11)

G4VPhysicalVolume* Em10DetectorConstruction::SetUpBarr90()
{
  fWorldSizeZ = 400.*cm;
  fWorldSizeR = 20.*cm;

  // Radiator and detector parameters

  fRadThickness = 0.019*mm;
  fGasGap       = 0.6*mm;
  foilGasRatio  = fRadThickness/(fRadThickness+fGasGap);

  fFoilNumber   = 350;

  fAbsorberThickness = 50.0*mm;

  fAbsorberRadius   = 100.*mm;
  fAbsorberZ        = 136.*cm;

  fWindowThick    = 51.0*micrometer;
  fElectrodeThick = 10.0*micrometer;
  fGapThick       =  10.0*cm;

  fDetThickness =  50.0*mm;
  fDetLength    = 200.0*cm;
  fDetGap       =   0.01*mm;

  fStartR       = 40*cm;
  fStartZ       = 100.0*mm;

  fModuleNumber = 1;

  // Preparation of mixed radiator material

  G4Material* CH2 = fMat->GetMaterial("CH2");
  G4Material* CO2 = fMat->GetMaterial("CO2");
  G4Material* Air   = fMat->GetMaterial("Air");
  G4Material* Al   = fMat->GetMaterial("Al");
  G4Material* Mylar = fMat->GetMaterial("Mylar");

  G4double foilDensity = 0.91*g/cm3;
  // CH21.39*g/cm3; // Mylar //  0.534*g/cm3; //Li
  G4double gasDensity  = 1.977*mg/cm3;
  // CO2 1.2928*mg/cm3; // Air //  0.178*mg/cm3; // He
 
  G4double totDensity  = foilDensity*foilGasRatio +
                                             gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;
 
  G4Material* radiatorMat = new G4Material("radiatorMat"  , totDensity,
                                                  2);
  radiatorMat->AddMaterial( CH2, fractionFoil );
  radiatorMat->AddMaterial( CO2, fractionGas  );

  // default materials of the detector and TR radiator

  fRadiatorMat =  radiatorMat;
  fFoilMat     = CH2;
  fGasMat      = CO2;  

  fWindowMat    = Mylar;
  fElectrodeMat = Al;

  fAbsorberMaterial = fMat->GetMaterial("Xe55He15CH4");
 
  fGapMat          = fAbsorberMaterial;

  fWorldMaterial    =  Air; // CO2; //

  fSolidWorld = new G4Box("World", fWorldSizeR,fWorldSizeR,fWorldSizeZ/2.);

  fLogicWorld = new G4LogicalVolume(fSolidWorld,  fWorldMaterial,  "World");

  fPhysicsWorld = new G4PVPlacement(0, G4ThreeVector(), "World",
                                 fLogicWorld, 0,  false, 0);

  // TR radiator envelope

  fRadThick = fFoilNumber*(fRadThickness + fGasGap) - fGasGap + fDetGap;

  fRadZ = fStartZ + 0.5*fRadThick;

  fSolidRadiator = new G4Box("Radiator",1.1*fAbsorberRadius ,
                              1.1*fAbsorberRadius,  0.5*fRadThick );

  fLogicRadiator = new G4LogicalVolume(fSolidRadiator, fRadiatorMat,
                                       "Radiator");

  fPhysicsRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,fRadZ),
                                     "Radiator", fLogicRadiator,
                                     fPhysicsWorld, false,        0 );

  // create region for window inside windowR for

  if( fRadRegion != 0 ) delete fRadRegion;
  if( fRadRegion == 0 ) fRadRegion = new G4Region("XTRradiator");
  fRadRegion->AddRootLogicalVolume(fLogicRadiator);

  fWindowZ = fStartZ + fRadThick + fWindowThick/2. + 15.0*mm;

  // G4Box* solidWindowR = new G4Box("WindowR",fAbsorberRadius+0.001,
  //                                         fAbsorberRadius+0.001,
  //                                         fWindowThick/2.+0.001  );

  // G4LogicalVolume* logicWindowR = new G4LogicalVolume(solidWindowR,
  //                                    fWorldMaterial, "WindowR");
  //
  //  G4VPhysicalVolume*    physiWindowR = new G4PVPlacement(0,
  //                       G4ThreeVector(0.,0.,fWindowZ),
  //                             "WindowR",logicWindowR,fPhysicsWorld,false,0);
  // window

  // G4Box* solidWindow = new G4Box("Window",fAbsorberRadius,
  //                                 fAbsorberRadius, fWindowThick/2.);

  // G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,
  //                                   fWindowMat, "Window");

  //  G4VPhysicalVolume*    physiWindow =
  //                        new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
  //                        "Window", logicWindow, physiWindowR, false, 0);

  fGapZ = fWindowZ + fWindowThick/2. + fGapThick/2. + 0.01*mm;

  fElectrodeZ = fGapZ + fGapThick/2. + fElectrodeThick/2. + 0.01*mm;

  // Absorber

  fAbsorberZ = fElectrodeZ + fElectrodeThick/2. +
                                             fAbsorberThickness/2. + 0.01*mm;

  fSolidAbsorber = new G4Box("Absorber", fAbsorberRadius,
                                 fAbsorberRadius, fAbsorberThickness/2.);

  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber, fAbsorberMaterial,
                                                "Absorber");
 
  fPhysicsAbsorber = new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsorberZ),
                                       "Absorber", fLogicAbsorber,
                                        fPhysicsWorld,  false,  0);

  if( fRegGasDet != 0 ) delete fRegGasDet;
  if( fRegGasDet == 0 ) fRegGasDet = new G4Region("XTRdEdxDetector");
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);

  // Sensitive Detectors: Absorber

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  PrintGeometryParameters();

  return fPhysicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::TestOld()
{
  //  G4double inch = 2.54*cm;
  // G4double  mil = inch/1000.0;
  //   G4double GetzstartAbs()           {return zstartAbs;};
  //  G4double GetzendAbs()             {return zendAbs;};
  // void ComputeCalorParameters();

  //  void SetGammaCut(G4double    cut){fGammaCut    = cut;};
  // void SetElectronCut(G4double cut){fElectronCut = cut;};
  //  void SetPositronCut(G4double cut){fPositronCut = cut;};
  // G4int fModelNumber ; // selection of parametrisation model1-10
  //   void SetAlphaPlate (G4double val){fAlphaPlate = val;};
  //   void SetAlphaGas   (G4double val){fAlphaGas   = val;};

  // G4double           fAlphaPlate ;
  // G4double           fAlphaGas ;

  // fAlphaPlate   = 160.0;
  // fAlphaGas     = 160.0;
  // fModelNumber  = 0;

  // create commands for interactive definition of the calorimeter

  // fGammaCut    = 23*mm;
  // fElectronCut = 23*mm; 
  // fPositronCut = 23*mm; 

  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //  G4int i, j ;
  // G4int j ;
  //  G4double zModule, zRadiator, rModule, rRadiator ;

  // complete the Calor parameters definition and Print

  //ComputeCalorParameters();

  // zRadiator ;
 
  // World

  // if(solidWorld) delete solidWorld ;
  // if(logicWorld) delete logicWorld ;
  // if(physiWorld) delete physiWorld ;
 
  //  if(solidRadiator) delete solidRadiator;
  //  if(logicRadiator) delete logicRadiator;
  //  if(physiRadiator) delete physiRadiator;

  //  radThick *= 1.02 ;

  //  if(fSolidRadSlice) delete fSolidRadSlice;
  //  if(fLogicRadSlice) delete fLogicRadSlice; 
  //  if(fPhysicRadSlice) delete fPhysicRadSlice;
  // fSolidRadSlice = new G4Box("RadSlice",fAbsorberRadius,
  //   fAbsorberRadius,0.5*fRadThickness );

  // fLogicRadSlice = new G4LogicalVolume(fSolidRadSlice,fRadiatorMat,
  //                                          "RadSlice",0,0,0);

//    for(j=0;j<fFoilNumber;j++)
//    {
//
//      zRadiator = zModule + j*(fRadThickness + fGasGap) ;
//      G4cout<<zRadiator/mm<<" mm"<<"\t" ;
//      //   G4cout<<"j = "<<j<<"\t" ;
// 
//      fPhysicRadSlice =
//          new G4PVPlacement(0,G4ThreeVector(0.,0.,zRadiator-zRad),
//                                         "RadSlice",fLogicRadSlice,
//                                          physiRadiator,false,j);
//     }
//  G4cout<<G4endl ;
 
    // fRadRegion->RemoveRootLogicalVolume(logicWindowR);
  // G4ProductionCuts* cutsR = 0;
    // cutsR = new G4ProductionCuts();
    // fRadRegion->SetProductionCuts(cutsR);

  // else  // Second time - get a cut object from region
  {
    // cutsR = fRadRegion->GetProductionCuts();
  }

  // cutsR->SetProductionCut(fGammaCut,"gamma");
  // cutsR->SetProductionCut(fElectronCut,"e-");
  // cutsR->SetProductionCut(fPositronCut,"e+");
  // G4Box* solidGap = new G4Box("Gap",fAbsorberRadius, fAbsorberRadius,
  //                                fGapThick/2.     ) ;
 
  // G4LogicalVolume* logicGap = new G4LogicalVolume(solidGap,fGapMat, "Gap");

  // G4VPhysicalVolume*    physiGap = new G4PVPlacement(0,
  //                                        G4ThreeVector(0.,0.,zGap),
  //                                    "Gap",logicGap,physiWorld,false,0);

  // G4Box* solidElectrode = new G4Box("Electrode",fAbsorberRadius,
  //                                  fAbsorberRadius, fElectrodeThick/2. );

  // G4LogicalVolume* logicElectrode = new G4LogicalVolume(solidElectrode,
  //                                     fElectrodeMat, "Electrode");

  //  G4VPhysicalVolume*    physiElectrode = new G4PVPlacement(0,
  //                                         G4ThreeVector(0.,0.,zElectrode),
  //                                    "Electrode",logicElectrode,
  //                                     physiWorld,false,0);
    //  if(solidAbsorber) delete solidAbsorber;
    //  if(logicAbsorber) delete logicAbsorber;
    //  if(physiAbsorber) delete physiAbsorber;
//   if (fAbsorberThickness > 0.)
//  {
//  }

    // fRegGasDet->RemoveRootLogicalVolume(logicAbsorber);
  // G4ProductionCuts* cuts = 0;
    // cuts = new G4ProductionCuts();
    //  fRegGasDet->SetProductionCuts(cuts);
  // else  // Second time - get a cut object from region
  {
    //  cuts = fRegGasDet->GetProductionCuts();
  }

  // cuts->SetProductionCut(fGammaCut,"gamma");
  // cuts->SetProductionCut(fElectronCut,"e-");
  // cuts->SetProductionCut(fPositronCut,"e+");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::PrintGeometryParameters()
{
  G4cout << "\n The  WORLD   is made of "
       << fWorldSizeZ/mm << "mm of " << fWorldMaterial->GetName();
  G4cout << ", the transverse size (R) of the world is " << 
                                         fWorldSizeR/mm << " mm. " << G4endl;
  G4cout << " The ABSORBER is made of "
       << fAbsorberThickness/mm << "mm of " << fAbsorberMaterial->GetName();
  G4cout << ", the transverse size (R) is " << fAbsorberRadius/mm << 
            " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " 
         << fAbsorberZ/mm << "  mm." << G4endl;

  G4cout<<"fRadZ = "<<fRadZ/mm<<" mm"<<G4endl;
 
  G4cout<<"fStartZ = "<<fStartZ/mm<<" mm"<<G4endl;

  G4cout<<"fRadThick = "<<fRadThick/mm<<" mm"<<G4endl;
  G4cout<<"fFoilNumber = "<<fFoilNumber<<G4endl;
  G4cout<<"fRadiatorMat = "<<fRadiatorMat->GetName()<<G4endl;
  G4cout<<"WorldMaterial = "<<fWorldMaterial->GetName()<<G4endl;
  //  G4cout<<"fAbsorberZ = "<<fAbsorberZ/mm<<" mm"<<G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name
  G4Material* pttoMaterial;

  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  {
    pttoMaterial = (*theMaterialTable)[J];
 
    if(pttoMaterial->GetName() == materialChoice)
    {
      fAbsorberMaterial = pttoMaterial;
      fLogicAbsorber->SetMaterial(pttoMaterial);
        // PrintCalorParameters();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetRadiatorMaterial(G4String materialChoice)
{
  // get the pointer to the material table

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name

  G4Material* pttoMaterial;
  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  {
    pttoMaterial = (*theMaterialTable)[J];

    if(pttoMaterial->GetName() == materialChoice)
    {
      fRadiatorMat = pttoMaterial;
//      fLogicRadSlice->SetMaterial(pttoMaterial);
      // PrintCalorParameters();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name
  G4Material* pttoMaterial;

  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  {
    pttoMaterial = (*theMaterialTable)[J];
 
    if(pttoMaterial->GetName() == materialChoice)
    {
      fWorldMaterial = pttoMaterial;
      fLogicWorld->SetMaterial(pttoMaterial);
       //  PrintCalorParameters();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  fAbsorberThickness = val;
  //  ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetRadiatorThickness(G4double val)
{
  // change XTR radiator thickness and recompute the calorimeter parameters
  fRadThickness = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetGasGapThickness(G4double val)
{
  // change XTR gas gap thickness and recompute the calorimeter parameters
  fGasGap = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberRadius(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fAbsorberRadius = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldChanged=true;
  fWorldSizeZ = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetWorldSizeR(G4double val)
{
  fWorldChanged=true;
  fWorldSizeR = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetAbsorberZpos(G4double val)
{
  fAbsorberZ  = val;
  // ComputeCalorParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::SetMagField(G4double)
{
  //apply a global uniform magnetic field along X axis

  /* *********************************************************

  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;             //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non null
  {
    magField = new G4UniformMagField(G4ThreeVector(fieldValue,0.,0.));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  }
  else
  {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }

  *************************************************************** */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetectorXTR());
}
