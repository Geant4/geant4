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
/// \file runAndEvent/RE01/src/RE01DetectorConstruction.cc
/// \brief Implementation of the RE01DetectorConstruction class
//
//

#include "RE01DetectorConstruction.hh"
#include "RE01TrackerSD.hh"
#include "RE01TrackerParametrisation.hh"
#include "RE01CalorimeterParametrisation.hh"
#include "RE01Field.hh"
#include "RE01RegionInformation.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01DetectorConstruction::RE01DetectorConstruction()
  : G4VUserDetectorConstruction()
{
#include "RE01DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01DetectorConstruction::~RE01DetectorConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* RE01DetectorConstruction::Construct()
{
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------

  G4double a, iz, density;
  G4String name, symbol;
  G4int nel;

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a);

  //  Material Information imported from NIST database.
  G4NistManager* NISTman = G4NistManager::Instance();
  G4Material* air     = NISTman->FindOrBuildMaterial("G4_AIR");
  G4Material* lead    = NISTman->FindOrBuildMaterial("G4_Pb");
  G4Material* arGas   = NISTman->FindOrBuildMaterial("G4_Ar");
  G4Material* silicon = NISTman->FindOrBuildMaterial("G4_Si");

  density = 1.032*g/cm3;
  G4Material* scinti = new G4Material(name="Scintillator", density, nel=2);
  scinti->AddElement(elC, 9);
  scinti->AddElement(elH, 10);

  //-------------------------------------------------------------------------
  // Detector geometry
  //-------------------------------------------------------------------------
  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",fExpHall_x,fExpHall_y,fExpHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,air,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log,"expHall_P",
                        0,false,0);
  G4VisAttributes* experimentalHallVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  experimentalHallVisAtt->SetForceWireframe(true);
  experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);
  G4Region* defaultRegion = (*(G4RegionStore::GetInstance()))[0];
  RE01RegionInformation* defaultRInfo = new RE01RegionInformation();
  defaultRInfo->SetWorld();
  defaultRInfo->Print();
  defaultRegion->SetUserInformation(defaultRInfo);

  //------------------------------ tracker
  G4VSolid * tracker_tubs
    = new G4Tubs("trkTubs_tubs",fTrkTubs_rmin,fTrkTubs_rmax,fTrkTubs_dz,
                 fTrkTubs_sphi,fTrkTubs_dphi);
  G4LogicalVolume * tracker_log
    = new G4LogicalVolume(tracker_tubs,arGas,"trackerT_L",0,0,0);
  // G4VPhysicalVolume * tracker_phys =
      new G4PVPlacement(0,G4ThreeVector(),tracker_log,"tracker_phys",
                        experimentalHall_log,false,0);
  G4VisAttributes* tracker_logVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  tracker_logVisAtt->SetForceWireframe(true);
  tracker_log->SetVisAttributes(tracker_logVisAtt);
  G4Region* trackerRegion = new G4Region("TrackerRegion");
  RE01RegionInformation* trackerInfo = new RE01RegionInformation();
  trackerInfo->SetTracker();
  trackerRegion->SetUserInformation(trackerInfo);
  tracker_log->SetRegion(trackerRegion);
  trackerRegion->AddRootLogicalVolume(tracker_log);

  //------------------------------ tracker layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * trackerLayer_tubs
    = new G4Tubs("trackerLayer_tubs",fTrkTubs_rmin,fTrkTubs_rmax,fTrkTubs_dz,
                 fTrkTubs_sphi,fTrkTubs_dphi);
  fTrackerLayer_log
    = new G4LogicalVolume(trackerLayer_tubs,silicon,"trackerB_L",0,0,0);
  G4VPVParameterisation * trackerParam
    = new RE01TrackerParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume *trackerLayer_phys =
      new G4PVParameterised("trackerLayer_phys",fTrackerLayer_log,tracker_log,
                           kXAxis, fNotrkLayers, trackerParam);
  G4VisAttributes* trackerLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  trackerLayer_logVisAtt->SetForceWireframe(true);
  fTrackerLayer_log->SetVisAttributes(trackerLayer_logVisAtt);

  //------------------------------ calorimeter
  G4VSolid * calorimeter_tubs
    = new G4Tubs("calorimeter_tubs",fCaloTubs_rmin,fCaloTubs_rmax,
                  fCaloTubs_dz,fCaloTubs_sphi,fCaloTubs_dphi);
  fCalorimeter_log
    = new G4LogicalVolume(calorimeter_tubs,scinti,"caloT_L",0,0,0);
  // G4VPhysicalVolume * calorimeter_phys =
  new G4PVPlacement(0,G4ThreeVector(),fCalorimeter_log,"caloM_P",
                        experimentalHall_log,false,0);
  G4VisAttributes* calorimeter_logVisATT
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  calorimeter_logVisATT->SetForceWireframe(true);
  fCalorimeter_log->SetVisAttributes(calorimeter_logVisATT);
  G4Region* calorimeterRegion = new G4Region("CalorimeterRegion");
  RE01RegionInformation* calorimeterInfo
    = new RE01RegionInformation();
  calorimeterInfo->SetCalorimeter();
  calorimeterRegion->SetUserInformation(calorimeterInfo);
  fCalorimeter_log->SetRegion(calorimeterRegion);
  calorimeterRegion->AddRootLogicalVolume(fCalorimeter_log);

  //------------------------------- Lead layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * caloLayer_tubs
    = new G4Tubs("caloLayer_tubs",fCaloRing_rmin,fCaloRing_rmax,
                  fCaloRing_dz,fCaloRing_sphi,fCaloRing_dphi);
  G4LogicalVolume * caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs,lead,"caloR_L",0,0,0);
  G4VPVParameterisation * calorimeterParam
    = new RE01CalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume * caloLayer_phys =
  new G4PVParameterised("caloLayer_phys",caloLayer_log,fCalorimeter_log,
                           kXAxis, fNocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_logVisAtt->SetForceWireframe(true);
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);


  return experimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01DetectorConstruction::ConstructSDandField()
{
  //------------------------------------------------------------------
  // Sensitive Detector
  //------------------------------------------------------------------
  
  G4String trackerSDname = "/mydet/tracker";
  RE01TrackerSD * trackerSD = new RE01TrackerSD(trackerSDname);
  G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD);
  SetSensitiveDetector(fTrackerLayer_log, trackerSD);

  // N.B. Calorimeter SD is defined in the parallel world.

  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------
  
  RE01Field* myField = new RE01Field;
  G4FieldManager* fieldMgr
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(myField);
  fieldMgr->CreateChordFinder(myField);
}

