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
/// \file analysis/A01/src/A01DetectorConstruction.cc
/// \brief Implementation of the A01DetectorConstruction class
//
// $Id$
// --------------------------------------------------------------
//

#include "A01DetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "A01DetectorConstMessenger.hh"
#include "A01MagneticField.hh"
#include "A01CellParameterisation.hh"
#include "A01Hodoscope.hh"
#include "A01DriftChamber.hh"
#include "A01EmCalorimeter.hh"

#include "G4PVReplica.hh"
#include "A01HadCalorimeter.hh"

A01DetectorConstruction::A01DetectorConstruction()
 : fAir(0), fArgonGas(0), fScintillator(0), fCsI(0), fLead(0),
   fWorldVisAtt(0), fMagneticVisAtt(0),
   fArmVisAtt(0), fHodoscopeVisAtt(0), fChamberVisAtt(0),
   fWirePlaneVisAtt(0), fEMcalorimeterVisAtt(0), fCellVisAtt(0),
   fHadCalorimeterVisAtt(0), fHadCalorimeterCellVisAtt(0),
   fArmAngle(30.*deg), fSecondArmPhys(0)

{
  fMessenger = new A01DetectorConstMessenger(this);
  fMagneticField = new A01MagneticField();
  fFieldMgr = new G4FieldManager();
  fArmRotation = new G4RotationMatrix();
  fArmRotation->rotateY(fArmAngle);
}

A01DetectorConstruction::~A01DetectorConstruction()
{
  delete fArmRotation;
  delete fMagneticField;
  delete fFieldMgr;
  delete fMessenger;

  DestroyMaterials();

  delete fWorldVisAtt;
  delete fMagneticVisAtt;
  delete fArmVisAtt;
  delete fHodoscopeVisAtt;
  delete fChamberVisAtt;
  delete fWirePlaneVisAtt;
  delete fEMcalorimeterVisAtt;
  delete fCellVisAtt;
  delete fHadCalorimeterVisAtt;
  delete fHadCalorimeterCellVisAtt;
}

G4VPhysicalVolume* A01DetectorConstruction::Construct()
{
  // All managed (deleted) by SDManager
  G4VSensitiveDetector* hodoscope1;
  G4VSensitiveDetector* hodoscope2;
  G4VSensitiveDetector* chamber1;
  G4VSensitiveDetector* chamber2;
  G4VSensitiveDetector* EMcalorimeter;
  G4VSensitiveDetector* HadCalorimeter;
  ConstructMaterials();

  // Local Magnetic Field
  static G4bool fieldIsInitialized = false;

  if(!fieldIsInitialized)
  {
    fFieldMgr->SetDetectorField(fMagneticField);
    fFieldMgr->CreateChordFinder(fMagneticField);
    fieldIsInitialized = true;
  }

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  G4VSolid* worldSolid = new G4Box("worldBox",10.*m,3.*m,10.*m);
  G4LogicalVolume* worldLogical
    = new G4LogicalVolume(worldSolid,fAir,"worldLogical",0,0,0);
  G4VPhysicalVolume* worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,0,0);

  // Tube with Local Magnetic field
   
  G4VSolid* magneticSolid = new G4Tubs("magneticTubs",0.,1.*m,1.*m,0.,360.*deg);
  G4NistManager* man = G4NistManager::Instance();
  G4Material* G4_Galactic = man->FindOrBuildMaterial("G4_Galactic");
  
   G4LogicalVolume* magneticLogical
    = new G4LogicalVolume(magneticSolid,G4_Galactic,"magneticLogical",fFieldMgr,0,0);
 
  // placement of Tube 
 
  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  fieldRot->rotateX(90.*deg);
  new G4PVPlacement(fieldRot,G4ThreeVector(),magneticLogical,
                    "magneticPhysical",worldLogical,0,0);

  // set "user limits" for drawing smooth curve
  G4UserLimits* userLimits = new G4UserLimits(5.0*cm);
  magneticLogical->SetUserLimits(userLimits);

  // first arm
  G4VSolid* firstArmSolid = new G4Box("firstArmBox",1.5*m,1.*m,3.*m);
  G4LogicalVolume* firstArmLogical
    = new G4LogicalVolume(firstArmSolid,fAir,"firstArmLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-5.*m),firstArmLogical,
                    "firstArmPhysical",worldLogical,0,0);

  // second arm
  G4VSolid* secondArmSolid = new G4Box("secondArmBox",2.*m,2.*m,3.5*m);
  G4LogicalVolume* secondArmLogical
    = new G4LogicalVolume(secondArmSolid,fAir,"secondArmLogical",0,0,0);
  G4double x = -5.*m * std::sin(fArmAngle);
  G4double z = 5.*m * std::cos(fArmAngle);
  fSecondArmPhys
    = new G4PVPlacement(fArmRotation,G4ThreeVector(x,0.,z),secondArmLogical,
                        "fSecondArmPhys",worldLogical,0,0);

  // hodoscopes in first arm
  G4VSolid* hodoscope1Solid = new G4Box("hodoscope1Box",5.*cm,20.*cm,0.5*cm);
  G4LogicalVolume* hodoscope1Logical
    = new G4LogicalVolume(hodoscope1Solid,fScintillator,"hodoscope1Logical",0,0,0);
  for(int i1=0;i1<15;i1++)
  {
    G4double x1 = (i1-7)*10.*cm;
    new G4PVPlacement(0,G4ThreeVector(x1,0.,-1.5*m),hodoscope1Logical,
                      "hodoscope1Physical",firstArmLogical,0,i1);
  }

  // drift chambers in first arm
  G4VSolid* chamber1Solid = new G4Box("chamber1Box",1.*m,30.*cm,1.*cm);
  G4LogicalVolume* chamber1Logical
    = new G4LogicalVolume(chamber1Solid,fArgonGas,"chamber1Logical",0,0,0);
  for(int j1=0;j1<5;j1++)
  {
    G4double z1 = (j1-2)*0.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),chamber1Logical,
                      "chamber1Physical",firstArmLogical,0,j1);
  }

  // "virtual" wire plane
  G4VSolid* wirePlane1Solid = new G4Box("wirePlane1Box",1.*m,30.*cm,0.1*mm);
  G4LogicalVolume* wirePlane1Logical
    = new G4LogicalVolume(wirePlane1Solid,fArgonGas,"wirePlane1Logical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),wirePlane1Logical,
                    "wirePlane1Physical",chamber1Logical,0,0);

  // hodoscopes in second arm
  G4VSolid* hodoscope2Solid = new G4Box("hodoscope2Box",5.*cm,20.*cm,0.5*cm);
  G4LogicalVolume* hodoscope2Logical
    = new G4LogicalVolume(hodoscope2Solid,fScintillator,"hodoscope2Logical",0,0,0);
  for(int i2=0;i2<25;i2++)
  {
    G4double x2 = (i2-12)*10.*cm;
    new G4PVPlacement(0,G4ThreeVector(x2,0.,0.),hodoscope2Logical,
                      "hodoscope2Physical",secondArmLogical,0,i2);
  }

  // drift chambers in second arm
  G4VSolid* chamber2Solid = new G4Box("chamber2Box",1.5*m,30.*cm,1.*cm);
  G4LogicalVolume* chamber2Logical
    = new G4LogicalVolume(chamber2Solid,fArgonGas,"chamber2Logical",0,0,0);
  for(int j2=0;j2<5;j2++)
  {
    G4double z2 = (j2-2)*0.5*m - 1.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z2),chamber2Logical,
                      "chamber2Physical",secondArmLogical,0,j2);
  }

  // "virtual" wire plane
  G4VSolid* wirePlane2Solid = new G4Box("wirePlane2Box",1.5*m,30.*cm,0.1*mm);
  G4LogicalVolume* wirePlane2Logical
    = new G4LogicalVolume(wirePlane2Solid,fArgonGas,"wirePlane2Logical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),wirePlane2Logical,
                    "wirePlane2Physical",chamber2Logical,0,0);

  // CsI calorimeter
  G4VSolid* EMcalorimeterSolid = new G4Box("EMcalorimeterBox",1.5*m,30.*cm,15.*cm);
  G4LogicalVolume* EMcalorimeterLogical
    = new G4LogicalVolume(EMcalorimeterSolid,fCsI,"EMcalorimeterLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*m),EMcalorimeterLogical,
                    "EMcalorimeterPhysical",secondArmLogical,0,0);

  // EMcalorimeter cells
  G4VSolid* cellSolid = new G4Box("cellBox",7.5*cm,7.5*cm,15.*cm);
  G4LogicalVolume* cellLogical
    = new G4LogicalVolume(cellSolid,fCsI,"cellLogical",0,0,0);
  G4VPVParameterisation* cellParam = new A01CellParameterisation();
  new G4PVParameterised("cellPhysical",cellLogical,EMcalorimeterLogical,
                         kXAxis,80,cellParam);

  // hadron calorimeter
  G4VSolid* HadCalorimeterSolid
    = new G4Box("HadCalorimeterBox",1.5*m,30.*cm,50.*cm);
  G4LogicalVolume* HadCalorimeterLogical
    = new G4LogicalVolume(HadCalorimeterSolid,fLead,"HadCalorimeterLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,3.*m),HadCalorimeterLogical,
                    "HadCalorimeterPhysical",secondArmLogical,0,0);

  // hadron calorimeter column
  G4VSolid* HadCalColumnSolid
    = new G4Box("HadCalColumnBox",15.*cm,30.*cm,50.*cm);
  G4LogicalVolume* HadCalColumnLogical
    = new G4LogicalVolume(HadCalColumnSolid,fLead,"HadCalColumnLogical",0,0,0);
  new G4PVReplica("HadCalColumnPhysical",HadCalColumnLogical,
                   HadCalorimeterLogical,kXAxis,10,30.*cm);

  // hadron calorimeter cell
  G4VSolid* HadCalCellSolid
    = new G4Box("HadCalCellBox",15.*cm,15.*cm,50.*cm);
  G4LogicalVolume* HadCalCellLogical
    = new G4LogicalVolume(HadCalCellSolid,fLead,"HadCalCellLogical",0,0,0);
  new G4PVReplica("HadCalCellPhysical",HadCalCellLogical,
                   HadCalColumnLogical,kYAxis,2,30.*cm);

  // hadron calorimeter layers
  G4VSolid* HadCalLayerSolid
    = new G4Box("HadCalLayerBox",15.*cm,15.*cm,2.5*cm);
  G4LogicalVolume* HadCalLayerLogical
    = new G4LogicalVolume(HadCalLayerSolid,fLead,"HadCalLayerLogical",0,0,0);
  new G4PVReplica("HadCalLayerPhysical",HadCalLayerLogical,
                  HadCalCellLogical,kZAxis,20,5.*cm);

  // scintillator plates
  G4VSolid* HadCalScintiSolid
    = new G4Box("HadCalScintiBox",15.*cm,15.*cm,0.5*cm);
  G4LogicalVolume* HadCalScintiLogical
    = new G4LogicalVolume(HadCalScintiSolid,fScintillator,"HadCalScintiLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*cm),HadCalScintiLogical,
                    "HadCalScintiPhysical",HadCalLayerLogical,0,0);

  // sensitive detectors -----------------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;

  hodoscope1 = new A01Hodoscope(SDname="/hodoscope1");
  SDman->AddNewDetector(hodoscope1);
  hodoscope1Logical->SetSensitiveDetector(hodoscope1);
  hodoscope2 = new A01Hodoscope(SDname="/hodoscope2");
  SDman->AddNewDetector(hodoscope2);
  hodoscope2Logical->SetSensitiveDetector(hodoscope2);

  chamber1 = new A01DriftChamber(SDname="/chamber1");
  SDman->AddNewDetector(chamber1);
  wirePlane1Logical->SetSensitiveDetector(chamber1);
  chamber2 = new A01DriftChamber(SDname="/chamber2");
  SDman->AddNewDetector(chamber2);
  wirePlane2Logical->SetSensitiveDetector(chamber2);

  EMcalorimeter = new A01EmCalorimeter(SDname="/EMcalorimeter");
  SDman->AddNewDetector(EMcalorimeter);
  cellLogical->SetSensitiveDetector(EMcalorimeter);

  HadCalorimeter = new A01HadCalorimeter(SDname="/HadCalorimeter");
  SDman->AddNewDetector(HadCalorimeter);
  HadCalScintiLogical->SetSensitiveDetector(HadCalorimeter);

  // visualization attributes ------------------------------------------------

  fWorldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fWorldVisAtt->SetVisibility(false);
  worldLogical->SetVisAttributes(fWorldVisAtt);

  fMagneticVisAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
  magneticLogical->SetVisAttributes(fMagneticVisAtt);

  fArmVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fArmVisAtt->SetVisibility(false);
  firstArmLogical->SetVisAttributes(fArmVisAtt);
  secondArmLogical->SetVisAttributes(fArmVisAtt);

  fHodoscopeVisAtt = new G4VisAttributes(G4Colour(0.8888,0.0,0.0));
  hodoscope1Logical->SetVisAttributes(fHodoscopeVisAtt);
  hodoscope2Logical->SetVisAttributes(fHodoscopeVisAtt);

  fChamberVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  chamber1Logical->SetVisAttributes(fChamberVisAtt);
  chamber2Logical->SetVisAttributes(fChamberVisAtt);

  fWirePlaneVisAtt = new G4VisAttributes(G4Colour(0.0,0.8888,0.0));
  fWirePlaneVisAtt->SetVisibility(false);
  wirePlane1Logical->SetVisAttributes(fWirePlaneVisAtt);
  wirePlane2Logical->SetVisAttributes(fWirePlaneVisAtt);

  fEMcalorimeterVisAtt = new G4VisAttributes(G4Colour(0.8888,0.8888,0.0));
  fEMcalorimeterVisAtt->SetVisibility(false);
  EMcalorimeterLogical->SetVisAttributes(fEMcalorimeterVisAtt);

  fCellVisAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
  cellLogical->SetVisAttributes(fCellVisAtt);

  fHadCalorimeterVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  HadCalorimeterLogical->SetVisAttributes(fHadCalorimeterVisAtt);
  fHadCalorimeterCellVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  fHadCalorimeterCellVisAtt->SetVisibility(false);
  HadCalColumnLogical->SetVisAttributes(fHadCalorimeterCellVisAtt);
  HadCalCellLogical->SetVisAttributes(fHadCalorimeterCellVisAtt);
  HadCalLayerLogical->SetVisAttributes(fHadCalorimeterCellVisAtt);
  HadCalScintiLogical->SetVisAttributes(fHadCalorimeterCellVisAtt);

  // return the world physical volume ----------------------------------------

  G4cout << G4endl << "The geometrical tree defined are : " << G4endl << G4endl;
  DumpGeometricalTree(worldPhysical);

  return worldPhysical;
}

void A01DetectorConstruction::ConstructMaterials()
{
  G4double a;
  G4double z;
  G4double density;
  G4double weightRatio;
  G4String name;
  G4String symbol;
  G4int nElem;

  // Argon gas
  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  fArgonGas = new G4Material(name="ArgonGas", z=18., a, density);

  // elements for mixtures and compounds
  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", z=8., a);
  a = 126.9*g/mole;
  G4Element* elI = new G4Element(name="Iodine", symbol="I", z=53., a);
  a = 132.9*g/mole;
  G4Element* elCs= new G4Element(name="Cesium", symbol="Cs", z=55., a);

  // Air
  density = 1.29*mg/cm3;
  fAir = new G4Material(name="Air", density, nElem=2);
  fAir->AddElement(elN, weightRatio=.7);
  fAir->AddElement(elO, weightRatio=.3);

  // Scintillator
  density = 1.032*g/cm3;
  fScintillator = new G4Material(name="Scintillator", density, nElem=2);
  fScintillator->AddElement(elC, 9);
  fScintillator->AddElement(elH, 10);

  // CsI
  density = 4.51*g/cm3;
  fCsI = new G4Material(name="CsI", density, nElem=2);
  fCsI->AddElement(elI, weightRatio=.5);
  fCsI->AddElement(elCs,weightRatio=.5);

  // Lead
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  fLead = new G4Material(name="Lead", z=82., a, density);

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void A01DetectorConstruction::DestroyMaterials()
{
  // Destroy all allocated elements and materials
  size_t i;
  G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();
  for(i=0;i<matTable->size();i++)
  { delete (*(matTable))[i]; }
  matTable->clear();
  G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable();
  for(i=0;i<elemTable->size();i++)
  { delete (*(elemTable))[i]; }
  elemTable->clear();
}

void A01DetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth)
{
  for(int isp=0;isp<depth;isp++)
  { G4cout << "  "; }
  G4cout << aVolume->GetName() << "[" << aVolume->GetCopyNo() << "] "
         << aVolume->GetLogicalVolume()->GetName() << " "
         << aVolume->GetLogicalVolume()->GetNoDaughters() << " "
         << aVolume->GetLogicalVolume()->GetMaterial()->GetName();
  if(aVolume->GetLogicalVolume()->GetSensitiveDetector())
  {
    G4cout << " " << aVolume->GetLogicalVolume()->GetSensitiveDetector()
                            ->GetFullPathName();
  }
  G4cout << G4endl;
  for(int i=0;i<aVolume->GetLogicalVolume()->GetNoDaughters();i++)
  { DumpGeometricalTree(aVolume->GetLogicalVolume()->GetDaughter(i),depth+1); }
}

void A01DetectorConstruction::SetArmAngle(G4double val)
{
  if(!fSecondArmPhys)
  {
    G4cerr << "Detector has not yet been constructed." << G4endl;
    return;
  }

  fArmAngle = val;
  *fArmRotation = G4RotationMatrix();  // make it unit vector
  fArmRotation->rotateY(fArmAngle);
  G4double x = -5.*m * std::sin(fArmAngle);
  G4double z = 5.*m * std::cos(fArmAngle);
  fSecondArmPhys->SetTranslation(G4ThreeVector(x,0.,z));

  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
