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
/// \file electromagnetic/TestEm8/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm8: Gaseous detector
//
// Created: 31.08.2010 V.Ivanchenko ob base of V.Grichine code
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
// 

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TargetSD.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4ProductionCuts.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

/////////////////////////////////////////////////////////////////////////////

DetectorConstruction::DetectorConstruction(PrimaryGeneratorAction* p)
  : fPhysWorld(0), fLogicWorld(0), fLogicWind(0), fLogicDet(0),
    fTargetSD(0), fRegGasDet(0), fPrimaryGenerator(p)
{
  fGasThickness = 23.0*mm;
  fGasRadius    = 10.*cm;

  fWindowThick  = 51.0*micrometer;

  DefineMaterials();

  fDetectorMessenger = new DetectorMessenger(this);

  G4double cut = 23.*mm;
  fGasDetectorCuts   = new G4ProductionCuts();
  fGasDetectorCuts->SetProductionCut(cut,"gamma");
  fGasDetectorCuts->SetProductionCut(cut,"e-");
  fGasDetectorCuts->SetProductionCut(cut,"e+");
  fGasDetectorCuts->SetProductionCut(cut,"proton");
}

//////////////////////////////////////////////////////////////////////////

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
  delete fGasDetectorCuts;
}

//////////////////////////////////////////////////////////////////////////////

void DetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials 
  G4String name, symbol ;          
  G4double density;  
  G4int nel; 
  G4int ncomponents; 
  G4double fractionmass;

  G4NistManager* manager = G4NistManager::Instance();
  //
  // define Elements
  //
  G4Element* elH  = manager->FindOrBuildElement(1);
  G4Element* elC  = manager->FindOrBuildElement(6);
  G4Element* elO  = manager->FindOrBuildElement(8);
  G4Element* elF  = manager->FindOrBuildElement(9);
  G4Element* elNe  = manager->FindOrBuildElement(10);
  G4Element* elXe = manager->FindOrBuildElement(54);
  //
  // simple gases at STP conditions 
  //
  G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  G4Material* Kr = manager->FindOrBuildMaterial("G4_Kr");
  G4Material* Xe     = manager->FindOrBuildMaterial("G4_Xe");
  // 
  // gases at STP conditions
  //
  G4Material* CarbonDioxide = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4Material* Mylar  = manager->FindOrBuildMaterial("G4_MYLAR");
  G4Material* Methane= manager->FindOrBuildMaterial("G4_METHANE");
  G4Material* Propane= manager->FindOrBuildMaterial("G4_PROPANE");
  G4Material* empty  = manager->FindOrBuildMaterial("G4_Galactic");

  // 93% Kr + 7% CH4, STP
  density = 3.491*mg/cm3 ;      
  G4Material* Kr7CH4 = new G4Material(name="Kr7CH4"  , density, 
                                      ncomponents=2);
  Kr7CH4->AddMaterial( Kr,       fractionmass = 0.986 ) ;
  Kr7CH4->AddMaterial( Methane,  fractionmass = 0.014 ) ;

  G4double TRT_Xe_density = 5.485*mg/cm3;
  G4Material* TRT_Xe = new G4Material(name="TRT_Xe", TRT_Xe_density, nel=1,
                                      kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_Xe->AddElement(elXe,1);

  G4double TRT_CO2_density = 1.842*mg/cm3;
  G4Material* TRT_CO2 = new G4Material(name="TRT_CO2", TRT_CO2_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CO2->AddElement(elC,1);
  TRT_CO2->AddElement(elO,2);

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CF4->AddElement(elC,1);
  TRT_CF4->AddElement(elF,4);

  // ATLAS TRT straw tube gas mixture (20 C, 1 atm)
  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 = new G4Material(name="XeCO2CF4", XeCO2CF4_density,
                                        ncomponents=3,
                                        kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);

  // C3H8,20 C, 2 atm
  density = 3.758*mg/cm3 ;
  G4Material* C3H8 = new G4Material(name="C3H8",density,nel=2) ;
  C3H8->AddElement(elC,3) ;
  C3H8->AddElement(elH,8) ;

  // 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1 atm 
  density = 4.9196*mg/cm3 ;
  G4Material* XeCH4C3H8 = new G4Material(name="XeCH4C3H8"  , 
                                  density,  ncomponents=3);
  XeCH4C3H8->AddMaterial( Xe,      fractionmass = 0.971 ) ;
  XeCH4C3H8->AddMaterial( Methane, fractionmass = 0.010 ) ;
  XeCH4C3H8->AddMaterial( Propane, fractionmass = 0.019 ) ;

  // 93% Ar + 7% CH4, STP
  density = 1.709*mg/cm3 ;      
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4", density, ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( Methane,  fractionmass = 0.029 ) ;

  // 80% Ar + 20% CO2, STP
  density = 1.8223*mg/cm3 ;      
  G4Material* Ar_80CO2_20 = new G4Material(name="ArCO2"  , density, 
                                           ncomponents=2);
  Ar_80CO2_20->AddMaterial( Argon,           fractionmass = 0.783 ) ;
  Ar_80CO2_20->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 80% Xe + 20% CO2, STP
  density = 5.0818*mg/cm3 ;      
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2", density, 
                                       ncomponents=2);
  Xe20CO2->AddMaterial( Xe,            fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide, fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP
  density = 3.601*mg/cm3 ;      
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2"  , density, 
                                       ncomponents=2);
  Kr20CO2->AddMaterial( Kr,            fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide, fractionmass = 0.11 ) ;

  // ALICE mixture TPC_Ne-CO2-2
  density = 0.939*mg/cm3 ;      
  G4Material* NeCO2 = new G4Material(name="TPC_Ne-CO2-2", density, 
                                            ncomponents=3);
  NeCO2->AddElement( elNe, fractionmass = 0.8039 ) ;
  NeCO2->AddElement( elO,  fractionmass = 0.1426 ) ;
  NeCO2->AddElement( elC,  fractionmass = 0.0535 ) ;
   
  fGasMat = XeCH4C3H8;
  fWindowMat = Mylar;
  fWorldMaterial = empty; 

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

/////////////////////////////////////////////////////////////////////////
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  if(fRegGasDet) { delete fRegGasDet; }
  fRegGasDet = new G4Region("GasDetector");
  fRegGasDet->SetProductionCuts(fGasDetectorCuts);

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4double contThick = fWindowThick*2 + fGasThickness;
  G4double contR     = fWindowThick*2 + fGasRadius;

  G4double worldSizeZ = contThick*1.2;
  G4double worldSizeR = contR*1.2;

  fPrimaryGenerator->SetPositionZ(-0.55*contThick);

  // Printout parameters
  G4cout << "\n The  WORLD   is made of " 
         << worldSizeZ/mm << "mm of " << fWorldMaterial->GetName() ;
  G4cout << ", the transverse size (R) of the world is " << worldSizeR/mm 
         << " mm. " << G4endl;
  G4cout << " The CONTAINER is made of " 
         << fWindowThick/mm << "mm of " << fWindowMat->GetName() << G4endl;
  G4cout << " The TARGET is made of " 
         << fGasThickness/mm << "mm of " << fGasMat->GetName() ;
  G4cout << ", the transverse size (R) is " << fGasRadius/mm << " mm. " << G4endl;
  G4cout << G4endl;
      
  // World
  G4Tubs* SolidWorld = new G4Tubs("World",                             
                                  0.,worldSizeR,worldSizeZ/2.,0.,CLHEP::twopi);
                   
  fLogicWorld = new G4LogicalVolume(SolidWorld, fWorldMaterial, "World");                
                                   
  fPhysWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(0.,0.,0.),     
                                 "World", 
                                 fLogicWorld,
                                 0,                        //its mother  volume
                                 false,                        //no boolean operation
                                 0);                        //copy number

  // Window
  G4Tubs* wind = new G4Tubs("Absorber",                
                            0.,contR,contThick/2.,0.,CLHEP::twopi); 

  fLogicWind = new G4LogicalVolume(wind, fWindowMat, "Window"); 

  G4PVPlacement* PhysWind = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "Window", 
                                              fLogicWind, fPhysWorld, false, 0);                
                                        
  // Detector volume
  G4Tubs* det = new G4Tubs("Gas", 0., fGasRadius, fGasThickness/2., 0., CLHEP::twopi);  

  fLogicDet = new G4LogicalVolume(det, fGasMat, "Gas"); 

  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "Gas", fLogicDet, PhysWind, false, 0);

  fRegGasDet->AddRootLogicalVolume(fLogicDet);

  // Sensitive Detectors:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(!fTargetSD)
    {
      fTargetSD = new TargetSD("GasSD");
      SDman->AddNewDetector( fTargetSD );
    }
  fLogicDet->SetSensitiveDetector(fTargetSD);

  // visualisation
  fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* color1 = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  fLogicWind->SetVisAttributes(color1);
  G4VisAttributes* color2 = new G4VisAttributes(G4Colour(0.0, 0.3, 0.7));
  fLogicDet->SetVisAttributes(color2);

  if(0.0 == fGasMat->GetIonisation()->GetMeanEnergyPerIonPair()) {
    SetPairEnergy(20*eV);
  }
  return fPhysWorld;
}

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetGasMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fGasMat) {
    G4cout << "### New target material: " << mat->GetName() << G4endl;
    fGasMat = mat;
    if(fLogicDet) { 
      fLogicDet->SetMaterial(mat); 
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetContainerMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fWindowMat) {
    G4cout << "### New material for container: " << mat->GetName() << G4endl;
    fWindowMat = mat;
    if(fLogicWind) { 
      fLogicWind->SetMaterial(mat); 
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetWorldMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fWorldMaterial) {
    G4cout << "### New World material: " << mat->GetName() << G4endl;
    fWorldMaterial = mat;
    if(fLogicWorld) { 
      fLogicWorld->SetMaterial(mat); 
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetGasThickness(G4double val)
{
  if(fGasThickness != val) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fGasThickness = val;
  }
}  

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetGasRadius(G4double val)
{
  if(fGasRadius != val) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fGasRadius = val;
  }
}  

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetContainerThickness(G4double val)
{
  if(fWindowThick != val) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fWindowThick = val;
  }
}  

///////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetPairEnergy(G4double val)
{
  if(val > 0.0) {
    fGasMat->GetIonisation()->SetMeanEnergyPerIonPair(val);
  }
}

////////////////////////////////////////////////////////////////////////////
