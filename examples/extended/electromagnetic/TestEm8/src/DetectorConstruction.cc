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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TargetSD.hh"

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
#include "TestParameters.hh"
#include "G4PionPlus.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fGasMat(nullptr), fWindowMat(nullptr), fWorldMaterial(nullptr),
    fSolidWorld(nullptr), fSolidContainer(nullptr), fSolidDetector(nullptr),
    fPhysWorld(nullptr), fLogicWorld(nullptr), fLogicContainer(nullptr), 
    fLogicDetector(nullptr),fRegGasDet(nullptr)
{
  fGasThickness = 23.0*mm;
  fGasRadius    = 10.*cm;
  fMaxStep      = DBL_MAX;

  fWindowThick  = 51.0*micrometer;

  DefineMaterials();

  fDetectorMessenger = new DetectorMessenger(this);

  G4double cut = 0.7*mm;
  fGasDetectorCuts   = new G4ProductionCuts();
  fGasDetectorCuts->SetProductionCut(cut,"gamma");
  fGasDetectorCuts->SetProductionCut(cut,"e-");
  fGasDetectorCuts->SetProductionCut(cut,"e+");
  fGasDetectorCuts->SetProductionCut(cut,"proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
  G4Element* elNe = manager->FindOrBuildElement(10);
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
  G4Material* CarbonDioxide = 
    manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4Material* Mylar  = manager->FindOrBuildMaterial("G4_MYLAR");
  G4Material* Methane= manager->FindOrBuildMaterial("G4_METHANE");
  G4Material* Propane= manager->FindOrBuildMaterial("G4_PROPANE");

  // propane at 10 atmospheres
  manager->ConstructNewGasMaterial("Propane10","G4_PROPANE",
                                   NTP_Temperature,10.*atmosphere);

  G4Material* empty  = manager->FindOrBuildMaterial("G4_Galactic");

  // 93% Kr + 7% CH4, STP
  density = 3.491*mg/cm3 ;
  G4Material* Kr7CH4 = 
    new G4Material(name="Kr7CH4"  , density, 
                   ncomponents=2);
  Kr7CH4->AddMaterial( Kr,       fractionmass = 0.986 ) ;
  Kr7CH4->AddMaterial( Methane,  fractionmass = 0.014 ) ;

  G4double TRT_Xe_density = 5.485*mg/cm3;
  G4Material* TRT_Xe = 
    new G4Material(name="TRT_Xe", TRT_Xe_density, nel=1,
                   kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_Xe->AddElement(elXe,1);

  G4double TRT_CO2_density = 1.842*mg/cm3;
  G4Material* TRT_CO2 = 
    new G4Material(name="TRT_CO2", TRT_CO2_density, nel=2,
                   kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CO2->AddElement(elC,1);
  TRT_CO2->AddElement(elO,2);

  // check alternative constructor
  std::vector<G4String> trtatom = {"C","O"};
  std::vector<G4int> trtnum  = {1, 2};
  manager->ConstructNewMaterial("TRT_CO2p",trtatom,trtnum,TRT_CO2_density,
                                true,kStateGas,NTP_Temperature,atmosphere);

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = 
    new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                   kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CF4->AddElement(elC,1);
  TRT_CF4->AddElement(elF,4);

  // ATLAS TRT straw tube gas mixture (20 C, 1 atm)
  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 = 
    new G4Material(name="XeCO2CF4", XeCO2CF4_density,
                   ncomponents=3,
                   kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);

  // C3H8,20 C, 2 atm
  density = 3.758*mg/cm3;
  G4Material* C3H8 = new G4Material(name="C3H8",density,nel=2,
                                    kStateGas,293.15*kelvin,2.*atmosphere);
  C3H8->AddElement(elC,3);
  C3H8->AddElement(elH,8);

  // The same material via different constructor
  std::vector<G4String> elmname = {"C","H"};
  std::vector<G4int>  atomnum = {3, 8};
  manager->ConstructNewIdealGasMaterial("C3H8p",elmname,atomnum,true,
                                        293.15*kelvin,2.*atmosphere);

  // 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1. atm 
  density = 4.9196*mg/cm3 ;
  G4Material* XeCH4C3H8 = 
    new G4Material(name="XeCH4C3H8"  , 
                   density,  ncomponents=3,
                   kStateGas,NTP_Temperature,1.*atmosphere);
  XeCH4C3H8->AddMaterial( Xe,      fractionmass = 0.971);
  XeCH4C3H8->AddMaterial( Methane, fractionmass = 0.010);
  XeCH4C3H8->AddMaterial( Propane, fractionmass = 0.019);

  // 93% Ar + 7% CH4, STP
  density = 1.709*mg/cm3 ;
  G4Material* Ar7CH4 = 
    new G4Material(name="Ar7CH4", density, ncomponents=2,
                   kStateGas,STP_Temperature,STP_Pressure);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( Methane,  fractionmass = 0.029 ) ;

  // 80% Ar + 20% CO2, STP
  density = 1.8223*mg/cm3 ;
  G4Material* Ar_80CO2_20 = 
    new G4Material(name="ArCO2"  , density, ncomponents=2,
                   kStateGas,STP_Temperature,STP_Pressure);
  Ar_80CO2_20->AddMaterial( Argon,           fractionmass = 0.783 ) ;
  Ar_80CO2_20->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 80% Xe + 20% CO2, STP
  density = 5.0818*mg/cm3 ;      
  G4Material* Xe20CO2 = 
    new G4Material(name="Xe20CO2", density, ncomponents=2,
                   kStateGas,STP_Temperature,STP_Pressure);
  Xe20CO2->AddMaterial( Xe,            fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide, fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP
  density = 3.601*mg/cm3 ;      
  G4Material* Kr20CO2 = 
    new G4Material(name="Kr20CO2", density, ncomponents=2,
                   kStateGas,STP_Temperature,STP_Pressure);
  Kr20CO2->AddMaterial( Kr,            fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide, fractionmass = 0.11 ) ;

  // ALICE mixture TPC_Ne-CO2-2
  density = 0.939*mg/cm3 ;      
  G4Material* NeCO2 = 
    new G4Material(name="TPC_Ne-CO2-2", density, ncomponents=3,
                   kStateGas,NTP_Temperature,1.*atmosphere);
  NeCO2->AddElement( elNe, fractionmass = 0.8039 ) ;
  NeCO2->AddElement( elO,  fractionmass = 0.1426 ) ;
  NeCO2->AddElement( elC,  fractionmass = 0.0535 ) ;

  // check alternative constructor
  std::vector<G4String> neatom = {"Ne","O","C"};
  std::vector<G4double> nefr   = {0.8039, 0.1426, 0.0536};
  manager->ConstructNewMaterial("TPC_Ne-CO2-2p",neatom,nefr,density,true,
                                kStateGas,NTP_Temperature,atmosphere);

  // ALICE TRD mixure 85% Xe + 15% CO2 NTP
  density = 4.9389*mg/cm3 ;      
  G4Material* Xe15CO2 = 
    new G4Material(name="Xe15CO2", density, ncomponents=2,
                   kStateGas,NTP_Temperature,atmosphere);
  Xe15CO2->AddMaterial( Xe,            fractionmass = 0.944 );
  Xe15CO2->AddMaterial( CarbonDioxide, fractionmass = 0.056 );
   
  fGasMat = XeCH4C3H8;
  fWindowMat = Mylar;
  fWorldMaterial = empty; 

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4double contThick = fWindowThick*2 + fGasThickness;
  G4double contR     = fWindowThick*2 + fGasRadius;

  G4double worldSizeZ = contThick*1.2;
  G4double worldSizeR = contR*1.2;

  TestParameters::GetPointer()->SetPositionZ(-0.55*contThick);

  // Printout parameters
  G4cout << "\n The  WORLD   is made of " 
         << worldSizeZ/mm << "mm of " << fWorldMaterial->GetName() ;
  G4cout << ", the transverse size (R) of the world is " << worldSizeR/mm 
         << " mm. " << G4endl;
  G4cout << " The CONTAINER is made of " 
         << fWindowThick/mm << "mm of " << fWindowMat->GetName() << G4endl;
  G4cout << " The TARGET is made of " 
         << fGasThickness/mm << "mm of " << fGasMat->GetName() ;
  G4cout << ", the transverse size (R) is " << fGasRadius/mm << " mm. " 
         << G4endl;
  G4cout << G4endl;
      
  // World
  fSolidWorld = 
    new G4Tubs("World",0.,worldSizeR,worldSizeZ/2.,0.,CLHEP::twopi);
                   
  fLogicWorld = new G4LogicalVolume(fSolidWorld, fWorldMaterial, "World");
                                   
  fPhysWorld = new G4PVPlacement(0,      
                                   G4ThreeVector(0.,0.,0.),     
                                 "World", 
                                 fLogicWorld,
                                 0,      
                                 false,  
                                 0);     

  // Window
  fSolidContainer = new G4Tubs("Absorber",                
                               0.,contR,contThick/2.,0.,CLHEP::twopi); 

  fLogicContainer = new G4LogicalVolume(fSolidContainer, fWindowMat, "Window"); 

  G4PVPlacement* PhysWind = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                                              "Window",  fLogicContainer,
                                              fPhysWorld, false, 0);
                                        
  // Detector volume
  fSolidDetector = new G4Tubs("Gas", 0., fGasRadius, fGasThickness/2.,
                              0., CLHEP::twopi); 

  fLogicDetector = new G4LogicalVolume(fSolidDetector, fGasMat, "Gas"); 

  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "Gas", fLogicDetector, 
                    PhysWind, false, 0);

  // defined gas detector region
  fRegGasDet = new G4Region("GasDetector");
  fRegGasDet->SetProductionCuts(fGasDetectorCuts);
  fRegGasDet->AddRootLogicalVolume(fLogicDetector);

  // visualisation
  fLogicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* color1 = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  fLogicContainer->SetVisAttributes(color1);
  G4VisAttributes* color2 = new G4VisAttributes(G4Colour(0.0, 0.3, 0.7));
  fLogicDetector->SetVisAttributes(color2);

  if(0.0 == fGasMat->GetIonisation()->GetMeanEnergyPerIonPair()) {
    SetPairEnergy(20*eV);
  }
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{ 
  auto sd = new TargetSD("GasSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  SetSensitiveDetector(fLogicDetector, sd); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fGasMat) {
    G4cout << "### New target material: " << mat->GetName() << G4endl;
    fGasMat = mat;
    if(fLogicDetector) { 
      fLogicDetector->SetMaterial(mat); 
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainerMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fWindowMat) {
    G4cout << "### New material for container: " << mat->GetName() << G4endl;
    fWindowMat = mat;
    if(fLogicContainer) { 
      fLogicContainer->SetMaterial(mat); 
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasThickness(G4double val)
{
  fGasThickness = val;
  if(fPhysWorld) { ChangeGeometry(); }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasRadius(G4double val)
{
  fGasRadius = val;
  if(fPhysWorld) { ChangeGeometry(); }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainerThickness(G4double val)
{
  fWindowThick = val;
  if(fPhysWorld) { ChangeGeometry(); }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPairEnergy(G4double val)
{
  if(val > 0.0) {
    fGasMat->GetIonisation()->SetMeanEnergyPerIonPair(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ChangeGeometry()
{
  G4double contThick = fWindowThick*2 + fGasThickness;
  G4double contR     = fWindowThick*2 + fGasRadius;

  G4double worldSizeZ = contThick*1.2;
  G4double worldSizeR = contR*1.2;

  TestParameters::GetPointer()->SetPositionZ(-0.55*contThick);

  fSolidWorld->SetOuterRadius(worldSizeR);
  fSolidWorld->SetZHalfLength(worldSizeZ*0.5);

  fSolidContainer->SetOuterRadius(contR);
  fSolidContainer->SetZHalfLength(contThick*0.5);

  fSolidDetector->SetOuterRadius(fGasRadius);
  fSolidDetector->SetZHalfLength(fGasThickness*0.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
