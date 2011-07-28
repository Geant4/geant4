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
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31DetectorConstruction -------
//
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31DetectorConstruction.hh"
#include "test31DetectorMessenger.hh"
#include "test31EventAction.hh"
#include "test31SD.hh"
#include "test31Histo.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4NistManager.hh"

#include "globals.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31DetectorConstruction::test31DetectorConstruction():
  AbsorberMaterial(0),
  WorldMaterial(0),
  solidWorld(0),
  logicWorld(0),
  physWorld(0),
  solidAbs(0),
  logicAbs(0),
  physAbs(0),
  magField(0),
  theEvent(0),
  myVerbose(0),
  detIsConstructed(false),
  nAbsSaved(0),
  nFirstEvtToDebug(-1),
  nLastEvtToDebug(-1)
{
  // Default parameter values of the calorimeter
  // corresponds to water test
  nameMatAbsorber   = G4String("G4_WATER");
  AbsorberThickness = 1.0*mm;    
  SizeXY            = 1000.0*mm;
  gap               = 0.0;
  NumberOfAbsorbers = 300;
  nameMatWorld      = G4String("G4_Galactic");
  nameMatGap        = G4String("G4_Galactic");
  WorldSizeZ        = 400.0*mm;
  maxDelta          = 10.0*MeV;

  DefineMaterials();
  test31Histo::GetPointer();

  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new test31DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31DetectorConstruction::~test31DetectorConstruction()
{ 
  delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* test31DetectorConstruction::Construct()
{
  return ConstructGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::DefineMaterials()
{ 
  if(myVerbose > 0) {
    G4cout << "test31DetectorConstruction: DefineMaterials starts" << G4endl;  
  } 

  G4int    ncomponents, natoms;
  G4Material* ma = 0;
  
  G4NistManager* man = G4NistManager::Instance();

  G4double density = 1.39*g/cm3;
  ma = new G4Material("Mylar"  , density, ncomponents=3);
  ma->AddElement(man->FindOrBuildElement("C"), natoms=10);
  ma->AddElement(man->FindOrBuildElement("H"), natoms=18);
  ma->AddElement(man->FindOrBuildElement("O"), natoms=5);
  if(ma) man->SetVerbose(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4VPhysicalVolume* test31DetectorConstruction::ConstructGeometry()
{
  if(myVerbose > 0) {
    G4cout << "test31DetectorConstruction: ConstructGeometry starts" << G4endl;
  } 

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  WorldMaterial = GetMaterial(nameMatWorld);
  AbsorberMaterial = GetMaterial(nameMatAbsorber);
  GapMaterial = GetMaterial(nameMatGap);
  ComputeGeomParameters();

  //     
  // World
  //
  solidWorld = new G4Box("World",SizeXY+1.0*mm,SizeXY+1.0*mm,WorldSizeZ);   
  logicWorld = new G4LogicalVolume(solidWorld,WorldMaterial,"World");
  physWorld  = new G4PVPlacement(0,G4ThreeVector(),"World",logicWorld,
                                 0,false,0);
  //                               
  // Absorber container
  // 
  G4double dz = 0.5*(AbsorberThickness + gap)*NumberOfAbsorbers;  
  G4Box* solid = new G4Box("Container",SizeXY,SizeXY,dz);
  G4LogicalVolume* lv = new G4LogicalVolume(solid,GapMaterial,"Container");
  G4VPhysicalVolume* pv = new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, dz),"Container",
					    lv,physWorld,false,0);
  //                               
  // Absorber
  // 
  solidAbs = new G4Box("Absorber",SizeXY,SizeXY,AbsorberThickness*0.5);
  logicAbs = new G4LogicalVolume(solidAbs,AbsorberMaterial,"Absorber");
  G4double z = 0.5*AbsorberThickness + gap - dz;

  for (G4int j=0; j<NumberOfAbsorbers; j++) {
  
    physAbs = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,z),
                                "Absorber",logicAbs,pv,false,j);
    z += (AbsorberThickness + gap); 
  }
  
  //                               
  // Sensitive Detectors: Absorber 
  //

  test31SD* calorimeterSD = new test31SD("test31");
  (G4SDManager::GetSDMpointer())->AddNewDetector( calorimeterSD );
  logicAbs->SetSensitiveDetector(calorimeterSD);

  PrintGeomParameters();  

  detIsConstructed = true;

  //
  //always return the physical World
  //
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::PrintGeomParameters()
{
  G4cout << "The  WORLD   is made of " 
         << " of " << WorldMaterial->GetName();
  G4cout << ". The transverse size (XY) of the world is " 
         << SizeXY/mm << " mm" << G4endl;
  G4cout << "The ABSORBER is made of " << NumberOfAbsorbers << " items of "
         << AbsorberThickness/mm  
         << " mm of " << AbsorberMaterial->GetName();
  G4cout << ". The transverse size (XY) is " 
         << SizeXY/mm << " mm" 
         << G4endl;
  G4cout << "The meanExc(eV)= " 
         << AbsorberMaterial->GetIonisation()->GetMeanExcitationEnergy()/eV
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Material* test31DetectorConstruction::GetMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(mat, false);
  if(!pttoMaterial) G4cout << "Find or BUild material " << mat << " fail " << G4endl;  
  if(detIsConstructed) MaterialIsChanged();
  return pttoMaterial;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::SetNumberOfAbsorbers(G4int val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  NumberOfAbsorbers = val;
  if(detIsConstructed) GeometryIsChanged();
}  
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
  if(detIsConstructed) GeometryIsChanged();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::SetAbsorberSizeXY(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  SizeXY = val;
  if(detIsConstructed) GeometryIsChanged();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::SetWorldSizeZ(G4double val)
{
  WorldSizeZ = val;
  if(detIsConstructed) GeometryIsChanged();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::SetGap(G4double val)
{
  gap = val;
  if(detIsConstructed) GeometryIsChanged();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::SetMagField(G4double fieldValue, G4int axis)
{
  // access to the field manager
  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if(magField) delete magField;		//delete the existing magn field
  
  // Create new field if >0
  if(fieldValue!=0.0) {

    G4ThreeVector B;
    // Choose direction of the field
    if(1 == axis) {
      B = G4ThreeVector(fieldValue,0.,0.);	
    } else if(2 == axis) {
      B = G4ThreeVector(0.,fieldValue,0.);	
    } else {
      B = G4ThreeVector(0.,0.,fieldValue);	
    }

    magField = new G4UniformMagField(B);        
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);

  // Set zero field
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::ComputeGeomParameters()
{
  // Compute derived parameters of the 1st absorber 
     
  if(WorldSizeZ < (AbsorberThickness + gap)*NumberOfAbsorbers)
     WorldSizeZ = (AbsorberThickness + gap)*NumberOfAbsorbers + 1.0*mm;

  test31Histo* h = test31Histo::GetPointer(); 
  h->SetNumberOfAbsorbers(NumberOfAbsorbers);
  h->SetAbsorberThickness(AbsorberThickness);
  h->SetNumAbsorbersSaved(nAbsSaved);
  h->SetGap(gap);
  h->SetAbsorberMaterial(AbsorberMaterial);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void test31DetectorConstruction::UpdateGeometry()
{
  (G4RunManager::GetRunManager())->DefineWorldVolume(ConstructGeometry());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::GeometryIsChanged()
{
  (G4RunManager::GetRunManager())->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31DetectorConstruction::MaterialIsChanged()
{
  (G4RunManager::GetRunManager())->CutOffHasBeenModified();
  (G4RunManager::GetRunManager())->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







