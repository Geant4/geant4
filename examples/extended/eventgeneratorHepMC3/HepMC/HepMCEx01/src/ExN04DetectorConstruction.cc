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
/// \file eventgenerator/HepMC/HepMCEx01/src/ExN04DetectorConstruction.cc
/// \brief Implementation of the ExN04DetectorConstruction class
//
//

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "ExN04CalorimeterParametrisation.hh"
#include "ExN04CalorimeterROGeometry.hh"
#include "ExN04CalorimeterSD.hh"
#include "ExN04DetectorConstruction.hh"
#include "ExN04Field.hh"
#include "ExN04MuonSD.hh"
#include "ExN04TrackerParametrisation.hh"
#include "ExN04TrackerSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04DetectorConstruction::ExN04DetectorConstruction()
 : G4VUserDetectorConstruction()
{
#include "ExN04DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04DetectorConstruction::~ExN04DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04DetectorConstruction::DefineMaterials()
{
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------
  G4NistManager* nistManager = G4NistManager::Instance();
  fAir = nistManager->FindOrBuildMaterial("G4_AIR");
  fLead = nistManager->FindOrBuildMaterial("G4_Pb");
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");

  G4double a, z, density;
  G4int nel;

  // Argon gas
  a= 39.95*g/mole;
  density= 1.782e-03*g/cm3;
  fAr= new G4Material("ArgonGas", z=18., a, density);

  // Scintillator
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  fScinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  fScinti-> AddElement(elC, 9);
  fScinti-> AddElement(elH, 10);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExN04DetectorConstruction::Construct()
{
  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------

  static G4bool fieldIsInitialized = false;
  if ( !fieldIsInitialized ) {
    ExN04Field* myField = new ExN04Field;
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()->
        GetFieldManager();
    fieldMgr-> SetDetectorField(myField);
    fieldMgr-> CreateChordFinder(myField);
    fieldIsInitialized = true;
  }

  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------
  DefineMaterials();

  //-------------------------------------------------------------------------
  // Detector geometry
  //-------------------------------------------------------------------------

  //------------------------------ experimental hall
  G4Box* experimentalHall_box =
         new G4Box("expHall_b", fexpHall_x, fexpHall_y, fexpHall_z);
  G4LogicalVolume* experimentalHall_log =
         new G4LogicalVolume(experimentalHall_box, fAir,"expHall_L", 0,0,0);
  G4VPhysicalVolume * experimentalHall_phys =
         new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log,
                           "expHall_P", 0, false,0);
  G4VisAttributes* experimentalHallVisAtt =
         new G4VisAttributes(G4Colour(1.,1.,1.));
  experimentalHallVisAtt-> SetForceWireframe(true);
  experimentalHall_log-> SetVisAttributes(experimentalHallVisAtt);

  //------------------------------ tracker
  G4VSolid* tracker_tubs
    = new G4Tubs("trkTubs_tubs", ftrkTubs_rmin, ftrkTubs_rmax, ftrkTubs_dz,
                 ftrkTubs_sphi, ftrkTubs_dphi);
  G4LogicalVolume* tracker_log
    = new G4LogicalVolume(tracker_tubs, fAr,"trackerT_L",0,0,0);
  // G4VPhysicalVolume * tracker_phys =
      new G4PVPlacement(0,G4ThreeVector(), tracker_log, "tracker_phys",
                        experimentalHall_log, false, 0);
  G4VisAttributes* tracker_logVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  tracker_logVisAtt->SetForceWireframe(true);
  tracker_log->SetVisAttributes(tracker_logVisAtt);

  //------------------------------ tracker layers
  // As an example for Parameterised volume
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid* trackerLayer_tubs
    = new G4Tubs("trackerLayer_tubs", ftrkTubs_rmin, ftrkTubs_rmax, ftrkTubs_dz,
                 ftrkTubs_sphi, ftrkTubs_dphi);
  G4LogicalVolume* trackerLayer_log
    = new G4LogicalVolume(trackerLayer_tubs, fSilicon,"trackerB_L",0,0,0);
  G4VPVParameterisation* trackerParam
    = new ExN04TrackerParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume *trackerLayer_phys =
      new G4PVParameterised("trackerLayer_phys", trackerLayer_log, tracker_log,
                           kXAxis, fnotrkLayers, trackerParam);
  G4VisAttributes* trackerLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  trackerLayer_logVisAtt->SetForceWireframe(true);
  trackerLayer_log->SetVisAttributes(trackerLayer_logVisAtt);

  //------------------------------ calorimeter
  G4VSolid* calorimeter_tubs
    = new G4Tubs("calorimeter_tubs", fcaloTubs_rmin, fcaloTubs_rmax,
                  fcaloTubs_dz, fcaloTubs_sphi, fcaloTubs_dphi);
  G4LogicalVolume* calorimeter_log
    = new G4LogicalVolume(calorimeter_tubs, fScinti, "caloT_L",0,0,0);
  // G4VPhysicalVolume * calorimeter_phys =
      new G4PVPlacement(0,G4ThreeVector(), calorimeter_log, "caloM_P",
                        experimentalHall_log, false,0);
  G4VisAttributes* calorimeter_logVisATT
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  calorimeter_logVisATT->SetForceWireframe(true);
  calorimeter_log->SetVisAttributes(calorimeter_logVisATT);

  //------------------------------- Lead layers
  // As an example for Parameterised volume
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid* caloLayer_tubs
    = new G4Tubs("caloLayer_tubs", fcaloRing_rmin, fcaloRing_rmax,
                  fcaloRing_dz, fcaloRing_sphi, fcaloRing_dphi);
  G4LogicalVolume* caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs, fLead, "caloR_L",0,0,0);
  G4VPVParameterisation* calorimeterParam
    = new ExN04CalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume * caloLayer_phys =
      new G4PVParameterised("caloLayer_phys",caloLayer_log,calorimeter_log,
                           kXAxis, fnocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_logVisAtt->SetForceWireframe(true);
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);

  //------------------------------ muon counters
  // As an example of CSG volumes with rotation
  G4VSolid* muoncounter_box
    = new G4Box("muoncounter_box", fmuBox_width, fmuBox_thick, fmuBox_length);
  G4LogicalVolume* muoncounter_log
    = new G4LogicalVolume(muoncounter_box, fScinti, "mucounter_L",0,0,0);
  G4VPhysicalVolume* muoncounter_phys;
  for( G4int i = 0; i < fnomucounter; i++ ) {
    G4double phi, x, y, z;
    phi = 360.*deg/fnomucounter*i;
    x = fmuBox_radius*std::sin(phi);
    y = fmuBox_radius*std::cos(phi);
    z = 0.*cm;
    G4RotationMatrix rm;
    rm.rotateZ(phi);
    muoncounter_phys
      = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(x,y,z)),
                          muoncounter_log, "muoncounter_P",
                          experimentalHall_log,false,i);
  }
  G4VisAttributes* muoncounter_logVisAtt
    = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  muoncounter_logVisAtt->SetForceWireframe(true);
  muoncounter_log->SetVisAttributes(muoncounter_logVisAtt);

  //------------------------------------------------------------------
  // Sensitive Detector
  //------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerSDname = "/mydet/tracker";
  ExN04TrackerSD * trackerSD = new ExN04TrackerSD(trackerSDname);
  SDman->AddNewDetector(trackerSD);
  trackerLayer_log->SetSensitiveDetector(trackerSD);

  G4String calorimeterSDname = "/mydet/calorimeter";
  ExN04CalorimeterSD* calorimeterSD = new ExN04CalorimeterSD(calorimeterSDname);
  G4String ROgeometryName = "CalorimeterROGeom";
  G4VReadOutGeometry* calRO = new ExN04CalorimeterROGeometry(ROgeometryName);
  calRO->BuildROGeometry();
  calorimeterSD->SetROgeometry(calRO);
  SDman->AddNewDetector(calorimeterSD);
  calorimeter_log->SetSensitiveDetector(calorimeterSD);

  G4String muonSDname = "/mydet/muon";
  ExN04MuonSD* muonSD = new ExN04MuonSD(muonSDname);
  SDman->AddNewDetector(muonSD);
  muoncounter_log->SetSensitiveDetector(muonSD);

  //------------------------------------------------------------------
  // Digitizer modules
  //------------------------------------------------------------------

  return experimentalHall_phys;
}
