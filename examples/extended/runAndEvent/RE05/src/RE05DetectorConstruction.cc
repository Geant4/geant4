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
/// \file RE05/src/RE05DetectorConstruction.cc
/// \brief Implementation of the RE05DetectorConstruction class
//

#include "RE05DetectorConstruction.hh"
#include "RE05TrackerSD.hh"
#include "RE05CalorimeterSD.hh"
#include "RE05MuonSD.hh"
#include "RE05TrackerParametrisation.hh"
#include "RE05CalorimeterParametrisation.hh"
#include "RE05Field.hh"

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05DetectorConstruction::RE05DetectorConstruction()
: G4VUserDetectorConstruction()
{

#include "RE05DetectorParameterDef.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05DetectorConstruction::~RE05DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05DetectorConstruction::DefineMaterials()
{
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------

  G4double a, z, density;
  G4int nel;

  G4Element* H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
  G4Element* C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen",   "O", z=8., a= 16.00*g/mole);

  fAir = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  fAir->AddElement(N, 70.*perCent);
  fAir->AddElement(O, 30.*perCent);

  fLead = 
  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  fAr = 
  new G4Material("ArgonGas",z=18., a= 39.95*g/mole, density=1.782*mg/cm3);

  fSilicon = 
  new G4Material("Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);

  fScinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  fScinti->AddElement(C, 9);
  fScinti->AddElement(H, 10);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RE05DetectorConstruction::Construct()
{
  DefineMaterials();

  //-------------------------------------------------------------------------
  // Detector geometry
  //-------------------------------------------------------------------------

  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",fExpHall_x,fExpHall_y,fExpHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,fAir,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log,"expHall_P",
                        0,false,0);
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  //------------------------------ tracker
  G4VSolid * fTracker_tubs
    = new G4Tubs("trkTubs_tubs",fTrkTubs_rmin,fTrkTubs_rmax,fTrkTubs_dz,
                 fTrkTubs_sphi,fTrkTubs_dphi);
  G4LogicalVolume * fTracker_log
    = new G4LogicalVolume(fTracker_tubs,fAr,"trackerT_L",0,0,0);
  // G4VPhysicalVolume * fTracker_phys =
      new G4PVPlacement(0,G4ThreeVector(),fTracker_log,"tracker_phys",
                        experimentalHall_log,false,0);
  G4VisAttributes* fTracker_logVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  fTracker_log->SetVisAttributes(fTracker_logVisAtt);

  //------------------------------ tracker layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * trackerLayer_tubs
    = new G4Tubs("trackerLayer_tubs",fTrkTubs_rmin,fTrkTubs_rmax,fTrkTubs_dz,
                 fTrkTubs_sphi,fTrkTubs_dphi);
  G4LogicalVolume * trackerLayer_log
    = new G4LogicalVolume(trackerLayer_tubs,fSilicon,"trackerB_L",0,0,0);
  G4VPVParameterisation * trackerParam
    = new RE05TrackerParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume *trackerLayer_phys =
      new G4PVParameterised("trackerLayer_phys",trackerLayer_log,fTracker_log,
                           kXAxis, fNotrkLayers, trackerParam);
  G4VisAttributes* trackerLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  trackerLayer_logVisAtt->SetForceWireframe(true);
  trackerLayer_log->SetVisAttributes(trackerLayer_logVisAtt);

  //------------------------------ calorimeter
  G4VSolid * calorimeter_tubs
    = new G4Tubs("calorimeter_tubs",fCaloTubs_rmin,fCaloTubs_rmax,
                  fCaloTubs_dz,fCaloTubs_sphi,fCaloTubs_dphi);
  G4LogicalVolume * calorimeter_log
    = new G4LogicalVolume(calorimeter_tubs,fScinti,"caloT_L",0,0,0);
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
    = new G4Tubs("caloLayer_tubs",fCaloRing_rmin,fCaloRing_rmax,
                  fCaloRing_dz,fCaloRing_sphi,fCaloRing_dphi);
  G4LogicalVolume * caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs,fLead,"caloR_L",0,0,0);
  G4VPVParameterisation * calorimeterParam
    = new RE05CalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume * caloLayer_phys =
      new G4PVParameterised("caloLayer_phys",caloLayer_log,calorimeter_log,
                           kXAxis, fNocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);

  //------------------------------ muon counters
  // As an example of CSG volumes with rotation
  G4VSolid * muoncounter_box
    = new G4Box("muoncounter_box",fMuBox_width,fMuBox_thick,
                fMuBox_length);
  G4LogicalVolume * muoncounter_log
    = new G4LogicalVolume(muoncounter_box,fScinti,"mucounter_L",0,0,0);
  for(int i=0; i<fNomucounter ; i++)
  {
    G4double phi, x, y, z;
    phi = 360.*deg/fNomucounter*i;
    x = fMuBox_radius*std::sin(phi);
    y = fMuBox_radius*std::cos(phi);
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

  return experimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05DetectorConstruction::ConstructSDandField()
{
  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------
  RE05Field* myField = new RE05Field;
  G4FieldManager* fieldMgr
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(myField);
  fieldMgr->CreateChordFinder(myField);

  //------------------------------------------------------------------
  // Sensitive Detectors
  //------------------------------------------------------------------
  G4String trackerSDname = "/mydet/tracker";
  RE05TrackerSD * trackerSD = new RE05TrackerSD(trackerSDname);
  G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD);
  SetSensitiveDetector("trackerB_L",trackerSD);

  G4String muonSDname = "/mydet/muon";
  RE05MuonSD * muonSD = new RE05MuonSD(muonSDname);
  G4SDManager::GetSDMpointer()->AddNewDetector(muonSD);
  SetSensitiveDetector("mucounter_L",muonSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
