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
// $Id: DetectorConstruction.cc,v 1.14 2010-12-21 19:43:08 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm9: Crystal calorimeter
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4TransportationManager.hh"

#include "G4GeometryManager.hh"
#include "G4FieldManager.hh"
#include "G4RunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4SDManager.hh"
#include "HcalSD.hh"
#include "HcalAbsSD.hh"
#include "EcalSD.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
  logicWorld(0),physiWorld(0),logicECal(0),logicCrystal(0)
{
  // create commands for interactive definition of the calorimeter
  detectorMessenger = new DetectorMessenger(this);
  magField = 0;
  
  //
  // define Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(1);
  G4Element* Zn  = manager->FindOrBuildElement(30);
  G4Element* Cu  = manager->FindOrBuildElement(29);
  //G4Material* Al  = manager->FindOrBuildMaterial("G4_Al");
  G4Material* Fe  = manager->FindOrBuildMaterial("G4_Fe");

  G4int natoms;
  G4double density, fractionmass;

  G4Material* CuZn = new G4Material("CuZn", density= 8.83*g/cm3, natoms=2);
  CuZn->AddElement(Cu, fractionmass = 0.7);
  CuZn->AddElement(Zn, fractionmass = 0.3);
  worldMaterial = manager->FindOrBuildMaterial("G4_AIR");
  ecalMaterial  = manager->FindOrBuildMaterial("G4_PbWO4");
  ecalMaterial->GetIonisation()->SetBirksConstant(.008415*mm/MeV);
  scinMaterial = manager->FindOrBuildMaterial("G4_POLYSTYRENE");

  // Geometry parameters of ECAl
  crystalLength = 230.0*mm;
  crystalWidth  = 22.0*mm;
  gap = .01*mm;

  // Geometry parameters of HCAl
  G4double absThick[17]={12.977, 7.815, 4.932, 4.932, 4.932, 4.932, 4.932, 4.932, 4.932, 4.932, 5.371, 5.516, 5.516, 5.516, 5.516, 5.516, 9.734};
  G4double scThick[17]={9, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 9}; 

  //  angularFactor = 1.4;
  angularFactor = 1.0;

  for(G4int i=0; i<17; i++){ 
    absorThickness[i]= absThick[i]*cm*angularFactor;
    scinThickness[i] = scThick[i]*mm*angularFactor;
    absorMaterial[i] = CuZn;
  }
  absorMaterial[0] = Fe;
  absorMaterial[1] = Fe;
  absorMaterial[16]= Fe;
  
  hcalWidth = 460.*mm;

  regionHCAL = 0;
  cutsHCAL = new G4ProductionCuts();
  cutsHCAL->SetProductionCut(1.0*mm);

  buildPreShower = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
  delete magField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()  
{
  // Compute derived geometry parameters

  ecalWidth = 7.*crystalWidth + 6.*gap;

  worldXY = std::max(1.5*hcalWidth, 1.2*ecalWidth);

  hcalThickness = 0.0;

  for(G4int i=0; i<17; i++){ 
    hcalThickness += (2*gap + absorThickness[i] + scinThickness[i]);
  }

  G4double coeff = 1.2;
  if(buildPreShower) coeff = 2.2;
  worldZ = coeff*(hcalThickness + crystalLength);
  posCenterHcalZ = crystalLength*.5; 
  posCenterEcalZ = -hcalThickness*.5;
  posCenterPreShowerZ = posCenterEcalZ - crystalLength;
  
  // Cleanup old geometry
  if(regionHCAL) {
    delete regionHCAL;
    regionHCAL = 0;
  }

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  G4Box* solidW = new G4Box("World",worldXY*.5,worldXY*.5,worldZ*.5);
  logicWorld = new G4LogicalVolume(solidW,worldMaterial,"World");
  physiWorld = new G4PVPlacement(0,G4ThreeVector(),"World",logicWorld,0,false,0);
  
  ConstructECAL();
  ConstructHCAL();
  if(buildPreShower) ConstructPreShower();

  HistoManager* man = HistoManager::GetPointer();
  man->SetWorldLength(worldZ);

  G4cout << "### New geometry is constructed HCAL width(cm)= " 
	 << hcalWidth/cm << G4endl;
  if(man->GetVerbose() > 0)
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  // always return world
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructECAL()
{
  // ECAL envelope
  G4Box* solidE = new G4Box("VolECal",ecalWidth*0.5,ecalWidth*0.5,crystalLength*0.5);
  logicECal = new G4LogicalVolume(solidE,ecalMaterial,"VolECal");
  G4VPhysicalVolume* physE = new G4PVPlacement(0,G4ThreeVector(0.,0.,posCenterEcalZ),
					       "VolECal",logicECal,physiWorld,false,0);

  G4double airWidth = ecalWidth - 2*crystalWidth;

  G4Box* solidGap = new G4Box("gap",airWidth*0.5,airWidth*0.5,crystalLength*0.5);
  G4LogicalVolume* logicGap = new G4LogicalVolume(solidGap,worldMaterial,"gap");
  G4VPhysicalVolume* physGap = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
						 "VolECal",logicGap,physE,false,0);

  // ECAL Sensitive detector
  EcalSD* sd = new EcalSD("Ecal");
  (G4SDManager::GetSDMpointer())->AddNewDetector(sd);


  // Crystals
  G4Box* solidCrystal = new G4Box("Crystal",crystalWidth*0.5,crystalWidth*0.5,
				  crystalLength*0.5);
  logicCrystal = new G4LogicalVolume( solidCrystal,ecalMaterial,"Crystal");

  G4double x0 = -(crystalWidth + gap)*2.0;
  G4double y  = x0;
  G4double x;
  G4int k = 0;
  G4VPhysicalVolume* pv;
  G4int i,j;

  for (i=0; i<5; i++) {
    x  = x0;
    for (j=0; j<5; j++) {

      pv = new G4PVPlacement(0,G4ThreeVector(x,y,0.),"Ecal",logicCrystal,
			     physGap,false,k);
      k++;
      x += crystalWidth + gap;
    }
    y += crystalWidth + gap;
  }

  //Color
  logicGap-> SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0, 0.5, .5));
  logicECal->SetVisAttributes(regWcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(.5, .5, .5));
  logicCrystal->SetVisAttributes(regCcolor);
  logicCrystal->SetSensitiveDetector(sd);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructHCAL()
{
  if(!regionHCAL) {
    regionHCAL = new G4Region("HcalRegion");
    regionHCAL->SetProductionCuts(cutsHCAL);
  }

  // HCAL envelope
  G4Box* solidHcal = new G4Box("Hcal",hcalWidth*0.5,hcalWidth*0.5,hcalThickness*0.5); 
  G4LogicalVolume* logicHcal = new G4LogicalVolume(solidHcal,worldMaterial,"Hcal");
  G4VPhysicalVolume* physHcal = new G4PVPlacement(0,G4ThreeVector(0.,0.,posCenterHcalZ),
						  "HCal",logicHcal,physiWorld,false,0);
  G4Box* box;
  G4Box* minibox; 
  G4LogicalVolume* logicV;
  G4LogicalVolume* logicMiniV; 
  G4PVPlacement* phys;

  //color
  logicHcal-> SetVisAttributes(G4VisAttributes::Invisible);

  // running position of a layer
  G4double z = -hcalThickness*0.5;

  //color absorber
  G4VisAttributes* BoxAbColor = new G4VisAttributes(G4Colour(0.5, 0, 0.5));
  G4VisAttributes* miniboxAbColor = new G4VisAttributes(G4Colour(1., 0, 1.));
  //color scint
  G4VisAttributes* BoxScintColor = new G4VisAttributes(G4Colour(0.5, 0.5, 0));
  G4VisAttributes* miniboxScintColor = new G4VisAttributes(G4Colour(1., 0., 0));

  // Sensitive detector
  HcalSD* scinsd = new HcalSD("Hcal");
  (G4SDManager::GetSDMpointer())->AddNewDetector(scinsd);
  HcalAbsSD* abssd = new HcalAbsSD("HcalAbs");
  (G4SDManager::GetSDMpointer())->AddNewDetector(abssd);
  
  for (G4int k=0; k<17; k++) {
    z += 0.5*absorThickness[k] + gap;
    
    // Absorber
    box = new G4Box("Ab",hcalWidth*0.6,hcalWidth*0.6,absorThickness[k]*0.5);
    logicV = new G4LogicalVolume(box,absorMaterial[k],"ab");
    phys = new G4PVPlacement(0,		               //no rotation
			     G4ThreeVector(0.,0.,z),   //its position
			     "ab",	               //its name
			     logicV,     	       //its logical volume	
			     physHcal,        	       //its mother
			     false,                    //no boulean operat
			     k);                       //copy number
    
    minibox = new G4Box("Ab1",hcalWidth*0.5,hcalWidth*0.5,absorThickness[k]*0.5);
    logicMiniV = new G4LogicalVolume(minibox,absorMaterial[k],"ab1");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),"ab1",logicMiniV,phys,false,k);

    z += 0.5*(absorThickness[k] + scinThickness[k]) + gap;

    logicV->SetVisAttributes(BoxAbColor);    
    logicMiniV->SetVisAttributes(miniboxAbColor);
    logicMiniV->SetSensitiveDetector(abssd);
    regionHCAL->AddRootLogicalVolume(logicMiniV);
    
    // Scint
    box = new G4Box("Sc",hcalWidth*0.6,hcalWidth*0.6,scinThickness[k]*0.5);
    logicV = new G4LogicalVolume(box,scinMaterial,"sc");
    phys = new G4PVPlacement(0,		               //no rotation
			     G4ThreeVector(0.,0.,z),   //its position
			     "sc",	               //its name
			     logicV,     	       //its logical volume	
			     physHcal,        	       //its mother
			     false,                    //no boulean operat
			     k);                       //copy number

    minibox = new G4Box("Sc1",hcalWidth*0.5,hcalWidth*0.5,scinThickness[k]*0.5);
    logicMiniV = new G4LogicalVolume(minibox,scinMaterial,"sc1");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),"sc1",logicMiniV,phys,false,k);

    logicV->SetVisAttributes(BoxScintColor);
    logicMiniV->SetVisAttributes(miniboxScintColor);
    logicMiniV->SetSensitiveDetector(scinsd);
    
    z += 0.5*scinThickness[k];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructPreShower()
{
  G4NistManager* manager = G4NistManager::Instance();
  const std::vector<G4String>& mat = manager->GetNistMaterialNames();
  G4int nmat  = mat.size();
  G4cout << "### DetectorConstruction::ConstructPreShower() added " 
	 << nmat << " materials " << G4endl;
  G4double dz = crystalLength/G4double(nmat);
  G4double z0 = posCenterPreShowerZ - crystalLength*0.5 - dz*0.5;
  G4Box* solid = new G4Box("PreShower",ecalWidth*0.5,ecalWidth*0.5,dz*0.5);
  G4LogicalVolume* logic;

  for(G4int i=0; i<nmat; i++) {

    z0 += dz;
    G4String matname = mat[i];
    G4Material* m = manager->FindOrBuildMaterial(matname);
    logic = new G4LogicalVolume(solid,m,matname);
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z0),matname,logic,physiWorld,false,0);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(mat, false);
  if (pttoMaterial && pttoMaterial != ecalMaterial) {
    ecalMaterial = pttoMaterial;
    if(logicCrystal) logicCrystal->SetMaterial(ecalMaterial);
    if(logicECal)    logicECal->SetMaterial(ecalMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalLength (G4double val)
{
  if(val > 0.0) {
    crystalLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCrystalWidth  (G4double val)
{
  if(val > 0.0) {
    crystalWidth = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapWidth  (G4double val)
{
  if(val > 0.0) {
    gap = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetHcalWidth  (G4double val)
{
  if(val > 0.0) {
    hcalWidth = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  // apply a global uniform magnetic field along X axis
  G4FieldManager* fieldMgr
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // create a new one if non nul
  if(fieldValue != 0.) {

    // delete the existing magn field
    if(magField) delete magField;	

    magField = new G4UniformMagField(G4ThreeVector(fieldValue,0.,0.));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);

    // zero field
  } else {
    delete magField;
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != worldMaterial) {
    worldMaterial = material;
    if(logicWorld) logicWorld->SetMaterial(worldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::SetBuildPreShower(G4bool val)
{
  buildPreShower = val;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
