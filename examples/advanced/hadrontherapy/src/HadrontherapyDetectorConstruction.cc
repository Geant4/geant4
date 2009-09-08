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
// HadrontherapyDetectorConstruction.cc
//
// See more at: http://workgroup.lngs.infn.it/geant4lns/
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova
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
#include "G4NistManager.hh"
#include "HadrontherapyDetectorROGeometry.hh"
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "PassiveProtonBeamLine.hh"
#include "HadrontherapyModulator.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction()
  : detectorSD(0), detectorROGeometry(0), 
    passiveProtonBeamLine(0), modulator(0),
    physicalTreatmentRoom(0),
    phantomPhysicalVolume(0), 
    detectorLogicalVolume(0), 
    detectorPhysicalVolume(0)
{
  // Messenger to change parameters of the geometry
  detectorMessenger = new HadrontherapyDetectorMessenger(this);

  // Detector sizes
  detectorSizeX = 20.*mm;
  detectorSizeY = 20.*mm;
  detectorSizeZ = 20.*mm;

  // Number of the detector voxels  
  numberOfVoxelsAlongX = 200;
  numberOfVoxelsAlongY = 1;
  numberOfVoxelsAlongZ = 1;
  
 }

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{ 
  if (detectorROGeometry) delete detectorROGeometry;  
  delete detectorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
//void HadrontherapyDetectorConstruction::ChangeTheBeamLine(const G4String& name)
//if (name == emName) return;
//if (name == "ProtonPassiveBeamLine")
//{
//ConstructPassiveProtonBeamLine();
//}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* HadrontherapyDetectorConstruction::Construct()
{ 

  ConstructPassiveProtonBeamLine();
  ConstructDetector();

  // Set the sensitive detector where the energy deposit is collected
  ConstructSensitiveDetector();
  return physicalTreatmentRoom;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::ConstructPassiveProtonBeamLine()
{ 
  // -----------------------------
  // Treatment room - World volume
  //------------------------------
  // Treatment room sizes
  const G4double worldX = 400.0 *cm;
  const G4double worldY = 400.0 *cm;
  const G4double worldZ = 400.0 *cm;
  G4bool isotopes = false;
 
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
  G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
                                                            airNist, 
                                                            "logicTreatmentRoom", 
							    0,0,0);
  physicalTreatmentRoom = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "physicalTreatmentRoom", 
					    logicTreatmentRoom, 
					    0,false,0);
 

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::Invisible);
 
  // Components of the Passive Proton Beam Line
  passiveProtonBeamLine = new PassiveProtonBeamLine(physicalTreatmentRoom);
  passiveProtonBeamLine -> HadrontherapyBeamLineSupport();
  passiveProtonBeamLine -> HadrontherapyBeamScatteringFoils();
  passiveProtonBeamLine -> HadrontherapyRangeShifter();
  passiveProtonBeamLine -> HadrontherapyBeamCollimators();
  passiveProtonBeamLine -> HadrontherapyBeamMonitoring();
  passiveProtonBeamLine -> HadrontherapyMOPIDetector();
  passiveProtonBeamLine -> HadrontherapyBeamNozzle();
  passiveProtonBeamLine -> HadrontherapyBeamFinalCollimator();

  // The following lines construc a typical modulator wheel inside the Passive Beam line.
  // Please remember to set the nodulator material (default is air, i.e. no modulator!) 
  // in the HadrontherapyModulator.cc file
  modulator = new HadrontherapyModulator();
  modulator -> BuildModulator(physicalTreatmentRoom);

  //----------------------------------------
  // Phantom:
  // A box used to approximate tissues
  //----------------------------------------

  G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
  G4Box* phantom = new G4Box("Phantom",20 *cm, 20 *cm, 20 *cm);
  G4LogicalVolume* phantomLogicalVolume = new G4LogicalVolume(phantom,	
							      waterNist, 
							      "phantomLog", 0, 0, 0);
  
  phantomPhysicalVolume = new G4PVPlacement(0,G4ThreeVector(200.*mm, 0.*mm, 0.*mm),
					    "phantomPhys",
					    phantomLogicalVolume,
					    physicalTreatmentRoom,
					    false,0);

  // Visualisation attributes of the patient
  red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
  red -> SetVisibility(true);
  red -> SetForceSolid(true);
  //red -> SetForceWireframe(true);
  phantomLogicalVolume -> SetVisAttributes(red);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::ConstructDetector()
{
  //-----------
  // Detector
  //-----------
  G4bool isotopes =  false; 
  G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
  G4Box* detector = new G4Box("Detector",detectorSizeX,detectorSizeY,detectorSizeZ);
  detectorLogicalVolume = new G4LogicalVolume(detector,
					      waterNist,
					      "DetectorLog",
					      0,0,0);
  
   G4double detectorXtranslation = -180.*mm;
  detectorPhysicalVolume = new G4PVPlacement(0,
					     G4ThreeVector(detectorXtranslation, 0.0 *mm, 0.0 *mm),
					     "DetectorPhys",
					     detectorLogicalVolume,
					     phantomPhysicalVolume,
					     false,0);
  
  // Visualisation attributes of the phantom
  skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
  skyBlue -> SetVisibility(true);
  skyBlue -> SetForceSolid(true);

  detectorLogicalVolume -> SetVisAttributes(skyBlue);
  
  // **************
  // Cut per Region
  // **************
  
  // A smaller cut is fixed in the phantom to calculate the energy deposit with the
  // required accuracy 
  G4Region* aRegion = new G4Region("DetectorLog");
  detectorLogicalVolume -> SetRegion(aRegion);
  aRegion -> AddRootLogicalVolume(detectorLogicalVolume);
}
/////////////////////////////////////////////////////////////////////////////
void  HadrontherapyDetectorConstruction::ConstructSensitiveDetector()
{  
  // Sensitive Detector and ReadOut geometry definition
  G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "Detector"; 

  if(!detectorSD)
    {
      // The sensitive detector is instantiated
      detectorSD = new HadrontherapyDetectorSD(sensitiveDetectorName);
      // The Read Out Geometry is instantiated
      G4String ROGeometryName = "DetectorROGeometry";
      detectorROGeometry = new HadrontherapyDetectorROGeometry(ROGeometryName,
							     detectorSizeX,
							     detectorSizeY,
							     detectorSizeZ,
							     numberOfVoxelsAlongX,
							     numberOfVoxelsAlongY,
							     numberOfVoxelsAlongZ);
      detectorROGeometry -> BuildROGeometry();
      detectorSD -> SetROgeometry(detectorROGeometry);
      sensitiveDetectorManager -> AddNewDetector(detectorSD);
      detectorLogicalVolume -> SetSensitiveDetector(detectorSD);
    }
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetModulatorAngle(G4double value)
{  
  modulator -> SetModulatorAngle(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetRangeShifterXPosition(G4double value)
{
  passiveProtonBeamLine -> SetRangeShifterXPosition(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetRSMaterial(G4String materialChoice)
{
  passiveProtonBeamLine -> SetRSMaterial(materialChoice);
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetRangeShifterXSize(G4double value)
{
  passiveProtonBeamLine -> SetRangeShifterXSize(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetFirstScatteringFoilSize(G4double value)
{
  passiveProtonBeamLine -> SetFirstScatteringFoilXSize(value);  
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetSecondScatteringFoilSize(G4double value)
{
  passiveProtonBeamLine -> SetSecondScatteringFoilXSize(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetOuterRadiusStopper(G4double value)
{
  passiveProtonBeamLine -> SetOuterRadiusStopper(value); 
  G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetInnerRadiusFinalCollimator(G4double value)
{
 passiveProtonBeamLine -> SetInnerRadiusFinalCollimator(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
