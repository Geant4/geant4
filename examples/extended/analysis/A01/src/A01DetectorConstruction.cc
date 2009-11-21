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
// $Id: A01DetectorConstruction.cc,v 1.10 2009-11-21 00:22:55 perl Exp $
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

#include "A01DetectorConstMessenger.hh"
#include "A01MagneticField.hh"
#include "A01CellParameterisation.hh"
#include "A01Hodoscope.hh"
#include "A01DriftChamber.hh"
#include "A01EmCalorimeter.hh"

#include "G4PVReplica.hh"
#include "A01HadCalorimeter.hh"

A01DetectorConstruction::A01DetectorConstruction()
 : air(0), argonGas(0), scintillator(0), CsI(0), lead(0),
   worldVisAtt(0), magneticVisAtt(0),
   armVisAtt(0), hodoscopeVisAtt(0), chamberVisAtt(0),
   wirePlaneVisAtt(0), EMcalorimeterVisAtt(0), cellVisAtt(0),
   HadCalorimeterVisAtt(0), HadCalorimeterCellVisAtt(0),
   armAngle(30.*deg), secondArmPhys(0)

{
  messenger = new A01DetectorConstMessenger(this);
  magneticField = new A01MagneticField();
  fieldMgr = new G4FieldManager();
  armRotation = new G4RotationMatrix();
  armRotation->rotateY(armAngle);
}

A01DetectorConstruction::~A01DetectorConstruction()
{
  delete armRotation;
  delete magneticField;
  delete fieldMgr;
  delete messenger;

  DestroyMaterials();

  delete worldVisAtt;
  delete magneticVisAtt;
  delete armVisAtt;
  delete hodoscopeVisAtt;
  delete chamberVisAtt;
  delete wirePlaneVisAtt;
  delete EMcalorimeterVisAtt;
  delete cellVisAtt;
  delete HadCalorimeterVisAtt;
  delete HadCalorimeterCellVisAtt;
}

G4VPhysicalVolume* A01DetectorConstruction::Construct()
{
  // All managed (deleted) by SDManager
  G4VSensitiveDetector* hodoscope1;
  G4VSensitiveDetector* hodoscope2;
  G4VSensitiveDetector* chamber1;
  G4VSensitiveDetector* chamber2;
  G4VSensitiveDetector* EMcalorimeter;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  G4VSensitiveDetector* HadCalorimeter;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ConstructMaterials();

  // Local Magnetic Field
  static G4bool fieldIsInitialized = false;

  if(!fieldIsInitialized)
  {
    fieldMgr->SetDetectorField(magneticField);
    fieldMgr->CreateChordFinder(magneticField);
    fieldIsInitialized = true;
  }

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  G4VSolid* worldSolid = new G4Box("worldBox",10.*m,3.*m,10.*m);
  G4LogicalVolume* worldLogical
    = new G4LogicalVolume(worldSolid,air,"worldLogical",0,0,0);
  G4VPhysicalVolume* worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,0,0);

  // Tube with Local Magnetic field
   
  G4VSolid* magneticSolid = new G4Tubs("magneticTubs",0.,1.*m,1.*m,0.,360.*deg);
  G4NistManager* man = G4NistManager::Instance();
  G4Material* G4_Galactic = man->FindOrBuildMaterial("G4_Galactic");
  
   G4LogicalVolume* magneticLogical
    = new G4LogicalVolume(magneticSolid,G4_Galactic,"magneticLogical",fieldMgr,0,0);
  //                                                                  ********
 
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
    = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-5.*m),firstArmLogical,
                    "firstArmPhysical",worldLogical,0,0);

  // second arm
  G4VSolid* secondArmSolid = new G4Box("secondArmBox",2.*m,2.*m,3.5*m);
  G4LogicalVolume* secondArmLogical
    = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical",0,0,0);
  G4double x = -5.*m * std::sin(armAngle);
  G4double z = 5.*m * std::cos(armAngle);
  secondArmPhys
    = new G4PVPlacement(armRotation,G4ThreeVector(x,0.,z),secondArmLogical,
                        "secondArmPhys",worldLogical,0,0);

  // hodoscopes in first arm
  G4VSolid* hodoscope1Solid = new G4Box("hodoscope1Box",5.*cm,20.*cm,0.5*cm);
  G4LogicalVolume* hodoscope1Logical
    = new G4LogicalVolume(hodoscope1Solid,scintillator,"hodoscope1Logical",0,0,0);
  for(int i1=0;i1<15;i1++)
  {
    G4double x1 = (i1-7)*10.*cm;
    new G4PVPlacement(0,G4ThreeVector(x1,0.,-1.5*m),hodoscope1Logical,
                      "hodoscope1Physical",firstArmLogical,0,i1);
  }

  // drift chambers in first arm
  G4VSolid* chamber1Solid = new G4Box("chamber1Box",1.*m,30.*cm,1.*cm);
  G4LogicalVolume* chamber1Logical
    = new G4LogicalVolume(chamber1Solid,argonGas,"chamber1Logical",0,0,0);
  for(int j1=0;j1<5;j1++)
  {
    G4double z1 = (j1-2)*0.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),chamber1Logical,
                      "chamber1Physical",firstArmLogical,0,j1);
  }

  // "virtual" wire plane
  G4VSolid* wirePlane1Solid = new G4Box("wirePlane1Box",1.*m,30.*cm,0.1*mm);
  G4LogicalVolume* wirePlane1Logical
    = new G4LogicalVolume(wirePlane1Solid,argonGas,"wirePlane1Logical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),wirePlane1Logical,
                    "wirePlane1Physical",chamber1Logical,0,0);

  // hodoscopes in second arm
  G4VSolid* hodoscope2Solid = new G4Box("hodoscope2Box",5.*cm,20.*cm,0.5*cm);
  G4LogicalVolume* hodoscope2Logical
    = new G4LogicalVolume(hodoscope2Solid,scintillator,"hodoscope2Logical",0,0,0);
  for(int i2=0;i2<25;i2++)
  {
    G4double x2 = (i2-12)*10.*cm;
    new G4PVPlacement(0,G4ThreeVector(x2,0.,0.),hodoscope2Logical,
                      "hodoscope2Physical",secondArmLogical,0,i2);
  }

  // drift chambers in second arm
  G4VSolid* chamber2Solid = new G4Box("chamber2Box",1.5*m,30.*cm,1.*cm);
  G4LogicalVolume* chamber2Logical
    = new G4LogicalVolume(chamber2Solid,argonGas,"chamber2Logical",0,0,0);
  for(int j2=0;j2<5;j2++)
  {
    G4double z2 = (j2-2)*0.5*m - 1.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z2),chamber2Logical,
                      "chamber2Physical",secondArmLogical,0,j2);
  }

  // "virtual" wire plane
  G4VSolid* wirePlane2Solid = new G4Box("wirePlane2Box",1.5*m,30.*cm,0.1*mm);
  G4LogicalVolume* wirePlane2Logical
    = new G4LogicalVolume(wirePlane2Solid,argonGas,"wirePlane2Logical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),wirePlane2Logical,
                    "wirePlane2Physical",chamber2Logical,0,0);

  // CsI calorimeter
  G4VSolid* EMcalorimeterSolid = new G4Box("EMcalorimeterBox",1.5*m,30.*cm,15.*cm);
  G4LogicalVolume* EMcalorimeterLogical
    = new G4LogicalVolume(EMcalorimeterSolid,CsI,"EMcalorimeterLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*m),EMcalorimeterLogical,
                    "EMcalorimeterPhysical",secondArmLogical,0,0);

  // EMcalorimeter cells
  G4VSolid* cellSolid = new G4Box("cellBox",7.5*cm,7.5*cm,15.*cm);
  G4LogicalVolume* cellLogical
    = new G4LogicalVolume(cellSolid,CsI,"cellLogical",0,0,0);
  G4VPVParameterisation* cellParam = new A01CellParameterisation();
  new G4PVParameterised("cellPhysical",cellLogical,EMcalorimeterLogical,
                         kXAxis,80,cellParam);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // hadron calorimeter
  G4VSolid* HadCalorimeterSolid
    = new G4Box("HadCalorimeterBox",1.5*m,30.*cm,50.*cm);
  G4LogicalVolume* HadCalorimeterLogical
    = new G4LogicalVolume(HadCalorimeterSolid,lead,"HadCalorimeterLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,3.*m),HadCalorimeterLogical,
                    "HadCalorimeterPhysical",secondArmLogical,0,0);

  // hadron calorimeter column
  G4VSolid* HadCalColumnSolid
    = new G4Box("HadCalColumnBox",15.*cm,30.*cm,50.*cm);
  G4LogicalVolume* HadCalColumnLogical
    = new G4LogicalVolume(HadCalColumnSolid,lead,"HadCalColumnLogical",0,0,0);
  new G4PVReplica("HadCalColumnPhysical",HadCalColumnLogical,
                   HadCalorimeterLogical,kXAxis,10,30.*cm);

  // hadron calorimeter cell
  G4VSolid* HadCalCellSolid
    = new G4Box("HadCalCellBox",15.*cm,15.*cm,50.*cm);
  G4LogicalVolume* HadCalCellLogical
    = new G4LogicalVolume(HadCalCellSolid,lead,"HadCalCellLogical",0,0,0);
  new G4PVReplica("HadCalCellPhysical",HadCalCellLogical,
                   HadCalColumnLogical,kYAxis,2,30.*cm);

  // hadron calorimeter layers
  G4VSolid* HadCalLayerSolid
    = new G4Box("HadCalLayerBox",15.*cm,15.*cm,2.5*cm);
  G4LogicalVolume* HadCalLayerLogical
    = new G4LogicalVolume(HadCalLayerSolid,lead,"HadCalLayerLogical",0,0,0);
  new G4PVReplica("HadCalLayerPhysical",HadCalLayerLogical,
                  HadCalCellLogical,kZAxis,20,5.*cm);

  // scintillator plates
  G4VSolid* HadCalScintiSolid
    = new G4Box("HadCalScintiBox",15.*cm,15.*cm,0.5*cm);
  G4LogicalVolume* HadCalScintiLogical
    = new G4LogicalVolume(HadCalScintiSolid,scintillator,"HadCalScintiLogical",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*cm),HadCalScintiLogical,
                    "HadCalScintiPhysical",HadCalLayerLogical,0,0);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

  worldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  worldVisAtt->SetVisibility(false);
  worldLogical->SetVisAttributes(worldVisAtt);

  magneticVisAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
  magneticLogical->SetVisAttributes(magneticVisAtt);

  armVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  armVisAtt->SetVisibility(false);
  firstArmLogical->SetVisAttributes(armVisAtt);
  secondArmLogical->SetVisAttributes(armVisAtt);

  hodoscopeVisAtt = new G4VisAttributes(G4Colour(0.8888,0.0,0.0));
  hodoscope1Logical->SetVisAttributes(hodoscopeVisAtt);
  hodoscope2Logical->SetVisAttributes(hodoscopeVisAtt);

  chamberVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  chamber1Logical->SetVisAttributes(chamberVisAtt);
  chamber2Logical->SetVisAttributes(chamberVisAtt);

  wirePlaneVisAtt = new G4VisAttributes(G4Colour(0.0,0.8888,0.0));
  wirePlaneVisAtt->SetVisibility(false);
  wirePlane1Logical->SetVisAttributes(wirePlaneVisAtt);
  wirePlane2Logical->SetVisAttributes(wirePlaneVisAtt);

  EMcalorimeterVisAtt = new G4VisAttributes(G4Colour(0.8888,0.8888,0.0));
  EMcalorimeterVisAtt->SetVisibility(false);
  EMcalorimeterLogical->SetVisAttributes(EMcalorimeterVisAtt);

  cellVisAtt = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
  cellLogical->SetVisAttributes(cellVisAtt);

  HadCalorimeterVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  HadCalorimeterLogical->SetVisAttributes(HadCalorimeterVisAtt);
  HadCalorimeterCellVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  HadCalorimeterCellVisAtt->SetVisibility(false);
  HadCalColumnLogical->SetVisAttributes(HadCalorimeterCellVisAtt);
  HadCalCellLogical->SetVisAttributes(HadCalorimeterCellVisAtt);
  HadCalLayerLogical->SetVisAttributes(HadCalorimeterCellVisAtt);
  HadCalScintiLogical->SetVisAttributes(HadCalorimeterCellVisAtt);

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
  argonGas = new G4Material(name="ArgonGas", z=18., a, density);

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
  air = new G4Material(name="Air", density, nElem=2);
  air->AddElement(elN, weightRatio=.7);
  air->AddElement(elO, weightRatio=.3);

  // Scintillator
  density = 1.032*g/cm3;
  scintillator = new G4Material(name="Scintillator", density, nElem=2);
  scintillator->AddElement(elC, 9);
  scintillator->AddElement(elH, 10);

  // CsI
  density = 4.51*g/cm3;
  CsI = new G4Material(name="CsI", density, nElem=2);
  CsI->AddElement(elI, weightRatio=.5);
  CsI->AddElement(elCs,weightRatio=.5);

  // Lead
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  lead = new G4Material(name="Lead", z=82., a, density);

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
  if(!secondArmPhys)
  {
    G4cerr << "Detector has not yet been constructed." << G4endl;
    return;
  }

  armAngle = val;
  *armRotation = G4RotationMatrix();  // make it unit vector
  armRotation->rotateY(armAngle);
  G4double x = -5.*m * std::sin(armAngle);
  G4double z = 5.*m * std::cos(armAngle);
  secondArmPhys->SetTranslation(G4ThreeVector(x,0.,z));

  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
