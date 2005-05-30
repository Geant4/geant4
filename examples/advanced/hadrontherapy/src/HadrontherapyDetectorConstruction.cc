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
// $Id: HadrontherapyDetectorConstruction.cc; Version 4.0 May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, G. Candiano, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyPhantomROGeometry.hh"
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4CSGSolid.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4TransportationManager.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyMaterial.hh"
#include "HadrontherapyBeamLine.hh"
#include "HadrontherapyModulator.hh"

HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction()
  : phantomSD(0), phantomROGeometry(0), beamLine(0),
    logicTreatmentRoom(0), physiTreatmentRoom(0), 
    PhantomLog(0), PhantomPhys(0), patientPhys(0)
{
// create commands for interactive definition of the calorimeter
  detectorMessenger = new HadrontherapyDetectorMessenger(this);

// PHANTOM SIZES
 phantomSizeX = 20.*mm;
 phantomSizeY = 20.*mm;
 phantomSizeZ = 20.*mm;

 numberOfVoxelsAlongX = 80;
 numberOfVoxelsAlongY = 80;
 numberOfVoxelsAlongZ = 80;

 ComputeDimVoxel();

 pMaterial = new HadrontherapyMaterial(); 
 modulator = new HadrontherapyModulator();
}

HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{ 
  delete modulator;
  delete pMaterial;

  if (phantomROGeometry) delete phantomROGeometry;
 delete detectorMessenger;
}

G4VPhysicalVolume* HadrontherapyDetectorConstruction::Construct()
{
 pMaterial -> DefineMaterials();

 ConstructBeamLine();
 ConstructPhantom();
 ConstructSensitiveDetector();
 
 return physiTreatmentRoom;
}

void HadrontherapyDetectorConstruction::ConstructBeamLine()
{
 G4Material* Air = pMaterial -> GetMat("Air") ;
 G4Material* Water = pMaterial -> GetMat("Water");

// ---------------------
// TREATMENT ROOM
//---------------------

 // TREATMENT ROOM SIZES
 G4double Worldx = 400.0 *cm;
 G4double Worldy = 400.0 *cm;
 G4double Worldz = 400.0 *cm;

 G4Box* TreatmentRoom = new G4Box("TreatmentRoom",Worldx,Worldy,Worldz);
 logicTreatmentRoom = new G4LogicalVolume(TreatmentRoom, Air, "logicTreatmentRoom", 0,0,0);
 physiTreatmentRoom = new G4PVPlacement(0,G4ThreeVector(),"physiTreatmentRoom", logicTreatmentRoom, 0,false,0);


 //Visualisation of the treatment room
 logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::Invisible);
 
 beamLine = new HadrontherapyBeamLine(physiTreatmentRoom);
 beamLine -> HadrontherapyBeamLineSupport();
 beamLine -> HadrontherapyBeamScatteringFoils();
 beamLine -> HadrontherapyBeamCollimators();
 beamLine -> HadrontherapyBeamMonitoring();
 beamLine -> HadrontherapyBeamNozzle();
 beamLine -> HadrontherapyBeamFinalCollimator();

 modulator -> BuildModulator(physiTreatmentRoom);

// Patient
 G4Box* patient = new G4Box("patient",20 *cm, 20 *cm, 20 *cm);
 G4LogicalVolume* patientLog = new G4LogicalVolume(patient,Water,
                                                   "patientLog",0,0,0);
 patientPhys = new G4PVPlacement(0,G4ThreeVector(0.0 *mm, 0.0 *mm, 0.0 *mm),
                                "patientPhys",
				 patientLog,physiTreatmentRoom,false,0); 

  G4VisAttributes * red_wire = new G4VisAttributes( G4Colour(1. ,0. ,0.));
  red_wire -> SetVisibility(true);
  red_wire -> SetForceWireframe(true);
  patientLog -> SetVisAttributes(red_wire); 
}

void HadrontherapyDetectorConstruction::ConstructPhantom()
{
  G4Colour  lblue   (0.0, 0.0, .75);

  G4Material* water = pMaterial->GetMat("Water");

  ComputeDimVoxel();

//----------------------
// WATER PHANTOM  

 G4Box* Phantom = new G4Box("Phantom",phantomSizeX,phantomSizeY,phantomSizeZ);

 PhantomLog = new G4LogicalVolume(Phantom,water,"PhantomLog",0,0,0);
 PhantomPhys = new G4PVPlacement(0,G4ThreeVector(-180.0 *mm, 0.0 *mm, 0.0 *mm),
                                 "PhantomPhys",
				  PhantomLog,patientPhys,false,0); 

 G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(lblue);
 simpleBoxVisAtt -> SetVisibility(true);
 simpleBoxVisAtt -> SetForceSolid(true);
 PhantomLog -> SetVisAttributes (simpleBoxVisAtt);

 // **************
 // Cut per Region
 // **************
 G4Region* aRegion = new G4Region("PhantomLog");
 PhantomLog -> SetRegion(aRegion);
 aRegion -> AddRootLogicalVolume(PhantomLog);
}

void  HadrontherapyDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "Phantom"; 

  if(!phantomSD)
  {
    phantomSD = new HadrontherapyPhantomSD(sensitiveDetectorName);
    G4String ROGeometryName = "PhantomROGeometry";
    phantomROGeometry = new HadrontherapyPhantomROGeometry(ROGeometryName,
							   phantomSizeX,
							   phantomSizeY,
							   phantomSizeZ,
							   numberOfVoxelsAlongX,
							   numberOfVoxelsAlongY,
							   numberOfVoxelsAlongZ);
    phantomROGeometry -> BuildROGeometry();
    phantomSD -> SetROgeometry(phantomROGeometry);
    pSDManager -> AddNewDetector(phantomSD);
    PhantomLog -> SetSensitiveDetector(phantomSD);
  }
}

void HadrontherapyDetectorConstruction::SetModulatorAngle(G4double val)
{  
  G4double modulatorAngle = val;
  G4cout << val/deg << G4endl;
  modulator -> SetModulatorAngle(modulatorAngle);
  G4RunManager::GetRunManager()-> GeometryHasBeenModified();
  //  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void HadrontherapyDetectorConstruction::SetRangeShifterXposition(G4double val)
{
  G4double value = val;
  beamLine -> setRangeShifterXPos(value);
  G4RunManager::GetRunManager()-> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetRangeShifterX(G4double val)
{
  G4double value = val;
  beamLine -> setRangeShifterX(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetFirstScatteringFoil(G4double val)
{
  G4double value = val;
  beamLine -> SetFirstScatteringFoil(value);  
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetSecondScatteringFoil(G4double val)
{
  G4double value = val;
  beamLine -> SetSecondScatteringFoil(value);
  G4RunManager::GetRunManager()-> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetOuterRadiusStopper(G4double val)
{
  G4double value = val;
  beamLine -> SetOuterRadiusStopper(value); 
  G4RunManager::GetRunManager()-> GeometryHasBeenModified(); 
 
  // Perche' si ridefinisce il sensitive detector?  
//ConstructSensitiveDetector();
  //G4RunManager::GetRunManager()-> GeometryHasBeenModified(); 
}

void HadrontherapyDetectorConstruction::SetInnerRadiusFinalCollimator(G4double val)
{
  G4double value = val;
  beamLine -> SetInnerRadiusFinalCollimator(value);
  G4RunManager::GetRunManager()-> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetRSMaterial(G4String materialChoice)
{
  G4String  newMaterial = materialChoice;
  beamLine -> SetRSMaterial(newMaterial);
}

/////////////////////////////////////////////////////////////////////////////
 void HadrontherapyDetectorConstruction::SetMaxStepSize(G4double val)
{
  // set the maximum allowed step size
  
 if (val <= DBL_MIN)
     { //G4cout << "\n --->warning from SetMaxStepSize: maxStep "
       // << val  << " out of range. Command refused" << G4endl;
      return;
     }
     userLimits->SetMaxAllowedStep(val);
     }

