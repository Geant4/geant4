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
// $Id: RE01DetectorConstruction.cc,v 1.3 2006-06-29 17:43:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


#include "RE01DetectorConstruction.hh"
#include "RE01TrackerSD.hh"
#include "RE01CalorimeterSD.hh"
#include "RE01CalorimeterROGeometry.hh"
#include "RE01TrackerParametrisation.hh"
#include "RE01CalorimeterParametrisation.hh"
#include "RE01Field.hh"
#include "RE01RegionInformation.hh"

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

RE01DetectorConstruction::RE01DetectorConstruction()
{

#include "RE01DetectorParameterDef.icc"

}

RE01DetectorConstruction::~RE01DetectorConstruction()
{;}

G4VPhysicalVolume* RE01DetectorConstruction::Construct()
{
  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------
/******************************************************************
  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)
  {
    RE01Field* myField = new RE01Field;
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
    fieldMgr->SetDetectorField(myField);
    fieldMgr->CreateChordFinder(myField);
    fieldIsInitialized = true;
  }
*******************************************************************/
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a);

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
  G4Material* Lead = new G4Material(name="Lead", z=82., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material * Silicon = new G4Material(name="Silicon", z=14., a, density);

  density = 1.032*g/cm3;
  G4Material* Scinti = new G4Material(name="Scintillator", density, nel=2);
  Scinti->AddElement(elC, 9);
  Scinti->AddElement(elH, 10);

  //-------------------------------------------------------------------------
  // Detector geometry
  //-------------------------------------------------------------------------

  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Air,"expHall_L",0,0,0);
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
  defaultRegion->SetUserInformation(defaultRInfo);

  //------------------------------ tracker
  G4VSolid * tracker_tubs
    = new G4Tubs("trkTubs_tubs",trkTubs_rmin,trkTubs_rmax,trkTubs_dz,
                 trkTubs_sphi,trkTubs_dphi);
  G4LogicalVolume * tracker_log
    = new G4LogicalVolume(tracker_tubs,Ar,"trackerT_L",0,0,0);
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
    = new G4Tubs("trackerLayer_tubs",trkTubs_rmin,trkTubs_rmax,trkTubs_dz,
                 trkTubs_sphi,trkTubs_dphi);
  G4LogicalVolume * trackerLayer_log
    = new G4LogicalVolume(trackerLayer_tubs,Silicon,"trackerB_L",0,0,0);
  G4VPVParameterisation * trackerParam
    = new RE01TrackerParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume *trackerLayer_phys =
      new G4PVParameterised("trackerLayer_phys",trackerLayer_log,tracker_log,
			   kXAxis, notrkLayers, trackerParam);
  G4VisAttributes* trackerLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  trackerLayer_logVisAtt->SetForceWireframe(true);
  trackerLayer_log->SetVisAttributes(trackerLayer_logVisAtt);

  //------------------------------ calorimeter
  G4VSolid * calorimeter_tubs
    = new G4Tubs("calorimeter_tubs",caloTubs_rmin,caloTubs_rmax,
		  caloTubs_dz,caloTubs_sphi,caloTubs_dphi);
  G4LogicalVolume * calorimeter_log
    = new G4LogicalVolume(calorimeter_tubs,Scinti,"caloT_L",0,0,0);
  // G4VPhysicalVolume * calorimeter_phys =
      new G4PVPlacement(0,G4ThreeVector(),calorimeter_log,"caloM_P",
			experimentalHall_log,false,0);
  G4VisAttributes* calorimeter_logVisATT
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  calorimeter_logVisATT->SetForceWireframe(true);
  calorimeter_log->SetVisAttributes(calorimeter_logVisATT);
  G4Region* calorimeterRegion = new G4Region("CalorimeterRegion");
  RE01RegionInformation* calorimeterInfo = new RE01RegionInformation();
  calorimeterInfo->SetCalorimeter();
  calorimeterRegion->SetUserInformation(calorimeterInfo);
  calorimeter_log->SetRegion(calorimeterRegion);
  calorimeterRegion->AddRootLogicalVolume(calorimeter_log);

  //------------------------------- Lead layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * caloLayer_tubs
    = new G4Tubs("caloLayer_tubs",caloRing_rmin,caloRing_rmax,
		  caloRing_dz,caloRing_sphi,caloRing_dphi);
  G4LogicalVolume * caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs,Lead,"caloR_L",0,0,0);
  G4VPVParameterisation * calorimeterParam
    = new RE01CalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume * caloLayer_phys =
      new G4PVParameterised("caloLayer_phys",caloLayer_log,calorimeter_log,
			   kXAxis, nocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_logVisAtt->SetForceWireframe(true);
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);

  //------------------------------------------------------------------
  // Sensitive Detector
  //------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerSDname = "/mydet/tracker";
  RE01TrackerSD * trackerSD = new RE01TrackerSD(trackerSDname);
  SDman->AddNewDetector(trackerSD);
  trackerLayer_log->SetSensitiveDetector(trackerSD);

  G4String calorimeterSDname = "/mydet/calorimeter";
  RE01CalorimeterSD * calorimeterSD = new RE01CalorimeterSD(calorimeterSDname);
  G4String ROgeometryName = "CalorimeterROGeom";
  G4VReadOutGeometry* calRO = new RE01CalorimeterROGeometry(ROgeometryName);
  calRO->BuildROGeometry();
  calorimeterSD->SetROgeometry(calRO);
  SDman->AddNewDetector(calorimeterSD);
  calorimeter_log->SetSensitiveDetector(calorimeterSD);

  return experimentalHall_phys;
}

