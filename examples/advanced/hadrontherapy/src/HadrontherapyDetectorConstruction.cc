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
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "HadrontherapyPhantomROGeometry.hh"
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyMaterial.hh"
#include "HadrontherapyBeamLine.hh"
#include "HadrontherapyModulator.hh"

HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction()
  : phantomSD(0), phantomROGeometry(0), beamLine(0), modulator(0),
    physicalTreatmentRoom(0),
    patientPhysicalVolume(0), 
    phantomLogicalVolume(0), 
    phantomPhysicalVolume(0)
{
  // Messenger to change parameters of the geometry
  detectorMessenger = new HadrontherapyDetectorMessenger(this);

  material = new HadrontherapyMaterial();

  // Phantom sizes
  phantomSizeX = 20.*mm;
  phantomSizeY = 20.*mm;
  phantomSizeZ = 20.*mm;

  // Number of the phantom voxels  
  numberOfVoxelsAlongX = 80;
  numberOfVoxelsAlongY = 80;
  numberOfVoxelsAlongZ = 80;
}

HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{ 
  delete material;
  if (phantomROGeometry) delete phantomROGeometry;  
  delete detectorMessenger;
}

G4VPhysicalVolume* HadrontherapyDetectorConstruction::Construct()
{  
  // Define the materials of the experimental set-up
  material -> DefineMaterials();
  
  // Define the geometry components
  ConstructBeamLine();
  ConstructPhantom();
  
  // Set the sensitive detector where the energy deposit is collected
  ConstructSensitiveDetector();
 
  return physicalTreatmentRoom;
}

void HadrontherapyDetectorConstruction::ConstructBeamLine()
{ 
  G4Material* air = material -> GetMat("Air") ;
  G4Material* water = material -> GetMat("Water");

  // ---------------------
  // Treatment room - World volume
  //---------------------

  // Treatment room sizes
  const G4double worldX = 400.0 *cm;
  const G4double worldY = 400.0 *cm;
  const G4double worldZ = 400.0 *cm;

  G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);

  G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
                                                            air, 
                                                            "logicTreatmentRoom", 
							    0,0,0);

  physicalTreatmentRoom = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "physicalTreatmentRoom", 
					    logicTreatmentRoom, 
					    0,false,0);


  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::Invisible);
 
  beamLine = new HadrontherapyBeamLine(physicalTreatmentRoom);
  beamLine -> HadrontherapyBeamLineSupport();
  beamLine -> HadrontherapyBeamScatteringFoils();
  beamLine -> HadrontherapyBeamCollimators();
  beamLine -> HadrontherapyBeamMonitoring();
  beamLine -> HadrontherapyBeamNozzle();
  beamLine -> HadrontherapyBeamFinalCollimator();

  modulator = new HadrontherapyModulator();
  modulator -> BuildModulator(physicalTreatmentRoom);

  // Patient - Mother volume of the phantom
  G4Box* patient = new G4Box("patient",20 *cm, 20 *cm, 20 *cm);

  G4LogicalVolume* patientLogicalVolume = new G4LogicalVolume(patient,
							      water, 
							      "patientLog", 0, 0, 0);

  patientPhysicalVolume = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),
					    "patientPhys",
					    patientLogicalVolume,
					    physicalTreatmentRoom,
					    false,0);
 
  // Visualisation attributes of the patient 
  G4VisAttributes * redWire = new G4VisAttributes(G4Colour(1. ,0. ,0.));
  redWire -> SetVisibility(true);
  redWire -> SetForceWireframe(true);
  patientLogicalVolume -> SetVisAttributes(redWire); 
}

void HadrontherapyDetectorConstruction::ConstructPhantom()
{
  G4Colour  lightBlue   (0.0, 0.0, .75);

  G4Material* water = material -> GetMat("Water");

  ComputeVoxelSize();

  //----------------------
  // Water phantom  
  //----------------------
  G4Box* phantom = new G4Box("Phantom",phantomSizeX,phantomSizeY,phantomSizeZ);

  phantomLogicalVolume = new G4LogicalVolume(phantom,
					     water,
					     "PhantomLog",
					     0,0,0);

  // Fixing the max step allowd in the phantom
  G4double maxStep = 0.02*cm;
  phantomLogicalVolume -> SetUserLimits(new G4UserLimits(maxStep));

  G4double phantomXtranslation = -180.*mm;
  phantomPhysicalVolume = new G4PVPlacement(0,
					    G4ThreeVector(phantomXtranslation, 0.0 *mm, 0.0 *mm),
					    "PhantomPhys",
					    phantomLogicalVolume,
					    patientPhysicalVolume,
					    false,0);
 
  // Visualisation attributes of the phantom
  G4VisAttributes* simpleBoxVisAttributes = new G4VisAttributes(lightBlue);
  simpleBoxVisAttributes -> SetVisibility(true);
  simpleBoxVisAttributes -> SetForceSolid(true);
  phantomLogicalVolume -> SetVisAttributes(simpleBoxVisAttributes);

  // **************
  // Cut per Region
  // **************
  
  // A smaller cut is fixed in the phantom to calculate the energy deposit with the
  // required accuracy 
  G4Region* aRegion = new G4Region("PhantomLog");
  phantomLogicalVolume -> SetRegion(aRegion);
  aRegion -> AddRootLogicalVolume(phantomLogicalVolume);
}

void  HadrontherapyDetectorConstruction::ConstructSensitiveDetector()
{  
  // Sensitive Detector and ReadOut geometry definition
  G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "Phantom"; 

  if(!phantomSD)
    {
      // The sensitive detector is instantiated
      phantomSD = new HadrontherapyPhantomSD(sensitiveDetectorName);
      
      // The Read Out Geometry is instantiated
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
      sensitiveDetectorManager -> AddNewDetector(phantomSD);
      phantomLogicalVolume -> SetSensitiveDetector(phantomSD);
    }
}

void HadrontherapyDetectorConstruction::SetModulatorAngle(G4double value)
{  
  modulator -> SetModulatorAngle(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetRangeShifterXPosition(G4double value)
{
  beamLine -> SetRangeShifterXPosition(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetRangeShifterXSize(G4double value)
{
  beamLine -> SetRangeShifterXSize(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetFirstScatteringFoilSize(G4double value)
{
  beamLine -> SetFirstScatteringFoilXSize(value);  
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetSecondScatteringFoilSize(G4double value)
{
  beamLine -> SetSecondScatteringFoilXSize(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetOuterRadiusStopper(G4double value)
{
  beamLine -> SetOuterRadiusStopper(value); 
  G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
}

void HadrontherapyDetectorConstruction::SetInnerRadiusFinalCollimator(G4double value)
{
  beamLine -> SetInnerRadiusFinalCollimator(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void HadrontherapyDetectorConstruction::SetRSMaterial(G4String materialChoice)
{
  beamLine -> SetRSMaterial(materialChoice);
}



