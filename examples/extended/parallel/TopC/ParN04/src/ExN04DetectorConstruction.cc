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
/// \file ExN04DetectorConstruction.cc
/// \brief Implementation of the ExN04DetectorConstruction class
//

#include "ExN04DetectorConstruction.hh"
#include "ExN04TrackerSD.hh"
#include "ExN04CalorimeterSD.hh"
#include "ExN04CalorimeterROGeometry.hh"
#include "ExN04MuonSD.hh"
#include "ExN04TrackerParametrisation.hh"
#include "ExN04CalorimeterParametrisation.hh"
#include "ExN04Field.hh"

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
#include "G4SystemOfUnits.hh"

ExN04DetectorConstruction::ExN04DetectorConstruction()
{

#include "ExN04DetectorParameterDef.icc"
  DefineMaterials();

}

ExN04DetectorConstruction::~ExN04DetectorConstruction()
{
  delete Scinti;
  delete Silicon;
  delete Ar;
  delete Lead;
  delete Air;

  delete O;
  delete N;
  delete C;
  delete H;
}

void ExN04DetectorConstruction::DefineMaterials()
{
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------

  G4double a, z, density;
  G4int nel;

  H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
  C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
  N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  O = new G4Element("Oxygen",   "O", z=8., a= 16.00*g/mole);

  Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  Lead = 
  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  Ar = 
  new G4Material("ArgonGas",z=18., a= 39.95*g/mole, density=1.782*mg/cm3);

  Silicon = 
  new G4Material("Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);

  Scinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  Scinti->AddElement(C, 9);
  Scinti->AddElement(H, 10);
}

G4VPhysicalVolume* ExN04DetectorConstruction::Construct()
{
  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------

  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)
  {
    ExN04Field* myField = new ExN04Field;
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
    fieldMgr->SetDetectorField(myField);
    fieldMgr->CreateChordFinder(myField);
    fieldIsInitialized = true;
  }


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
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

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
  tracker_log->SetVisAttributes(tracker_logVisAtt);

  //------------------------------ tracker layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * trackerLayer_tubs
    = new G4Tubs("trackerLayer_tubs",trkTubs_rmin,trkTubs_rmax,trkTubs_dz,
                 trkTubs_sphi,trkTubs_dphi);
  G4LogicalVolume * trackerLayer_log
    = new G4LogicalVolume(trackerLayer_tubs,Silicon,"trackerB_L",0,0,0);
  G4VPVParameterisation * trackerParam
    = new ExN04TrackerParametrisation;
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

  //------------------------------- Lead layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * caloLayer_tubs
    = new G4Tubs("caloLayer_tubs",caloRing_rmin,caloRing_rmax,
                  caloRing_dz,caloRing_sphi,caloRing_dphi);
  G4LogicalVolume * caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs,Lead,"caloR_L",0,0,0);
  G4VPVParameterisation * calorimeterParam
    = new ExN04CalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume * caloLayer_phys =
      new G4PVParameterised("caloLayer_phys",caloLayer_log,calorimeter_log,
                           kXAxis, nocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);

  //------------------------------ muon counters
  // As an example of CSG volumes with rotation
  G4VSolid * muoncounter_box
    = new G4Box("muoncounter_box",muBox_width,muBox_thick,
                muBox_length);
  G4LogicalVolume * muoncounter_log
    = new G4LogicalVolume(muoncounter_box,Scinti,"mucounter_L",0,0,0);
  for(int i=0; i<nomucounter ; i++)
  {
    G4double phi, x, y, z;
    phi = 360.*deg/nomucounter*i;
    x = muBox_radius*std::sin(phi);
    y = muBox_radius*std::cos(phi);
    z = 0.*cm;
    G4RotationMatrix rm;
    rm.rotateZ(phi);
    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(x,y,z)),
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
  ExN04CalorimeterSD * calorimeterSD = new ExN04CalorimeterSD(calorimeterSDname);
  G4String ROgeometryName = "CalorimeterROGeom";
  G4VReadOutGeometry* calRO = new ExN04CalorimeterROGeometry(ROgeometryName);
  calRO->BuildROGeometry();
  calRO->SetName(ROgeometryName);
  calorimeterSD->SetROgeometry(calRO);
  SDman->AddNewDetector(calorimeterSD);
  calorimeter_log->SetSensitiveDetector(calorimeterSD);

  G4String muonSDname = "/mydet/muon";
  ExN04MuonSD * muonSD = new ExN04MuonSD(muonSDname);
  SDman->AddNewDetector(muonSD);
  muoncounter_log->SetSensitiveDetector(muonSD);

  //------------------------------------------------------------------
  // Digitizer modules
  //------------------------------------------------------------------

  return experimentalHall_phys;
}

