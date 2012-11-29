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
/// \file analysis/N03Con/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fAbsorberMaterial(NULL),fGapMaterial(NULL),fDefaultMaterial(NULL),
 /*solidWorld(0),logicWorld(0),*/fPhysiWorld(NULL),
 /*solidCalor(0),logicCalor(0),physiCalor(0),*/
 /*solidLayer(0),logicLayer(0),physiLayer(0),*/
 /*solidAbsorber(0),logicAbsorber(0),*/fPhysiAbsorber(NULL),
 /*solidGap (0),logicGap (0),*/fPhysiGap (NULL),
 fMagField(NULL)
{
  // default parameter values of the calorimeter
  fAbsorberThickness = 10.*mm;
  fGapThickness      =  5.*mm;
  fNbOfLayers        = 10;
  fCalorSizeYZ       = 10.*cm;
  ComputeCalorParameters();
  
  // materials
  DefineMaterials();
  SetAbsorberMaterial("Lead");
  SetGapMaterial("liquidArgon");
  
  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  
G4int iz, n;                 //iz=number of protons  in an isotope; 
                             // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double abundance, fractionmass;

//
// define Elements
//

G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);

//
// define an Element from isotopes, by relative abundance 
//

G4Isotope* U5 = new G4Isotope("U235", iz=92, n=235, a=235.01*g/mole);
G4Isotope* U8 = new G4Isotope("U238", iz=92, n=238, a=238.03*g/mole);

G4Element* U  = new G4Element("enriched Uranium",symbol="U",ncomponents=2);
U->AddIsotope(U5, abundance= 90.*perCent);
U->AddIsotope(U8, abundance= 10.*perCent);

//
// define simple materials
//

new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

//
// define a material from elements.   case 1: chemical molecule
//

G4Material* H2O = 
new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);
// overwrite computed meanExcitationEnergy with ICRU recommended value 
H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

G4Material* Sci = 
new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
Sci->AddElement(C, natoms=9);
Sci->AddElement(H, natoms=10);

G4Material* Myl = 
new G4Material("Mylar", density= 1.397*g/cm3, ncomponents=3);
Myl->AddElement(C, natoms=10);
Myl->AddElement(H, natoms= 8);
Myl->AddElement(O, natoms= 4);

G4Material* SiO2 = 
new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
SiO2->AddElement(Si, natoms=1);
SiO2->AddElement(O , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

G4Material* Air = 
new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

G4Material* Aerog = 
new G4Material("Aerogel", density= 0.200*g/cm3, ncomponents=3);
Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
Aerog->AddElement (C   , fractionmass= 0.1*perCent);

//
// examples of gas in non STP conditions
//

G4Material* CO2 = 
new G4Material("CarbonicGas", density= 1.842*mg/cm3, ncomponents=2,
                              kStateGas, 325.*kelvin, 50.*atmosphere);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);
 
G4Material* steam = 
new G4Material("WaterSteam", density= 0.3*mg/cm3, ncomponents=1,
                             kStateGas, 500.*kelvin, 2.*atmosphere);
steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

G4Material* vacuum =
new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

G4Material* beam = 
new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                       kStateGas, STP_Temperature, 2.e-2*bar);
beam->AddMaterial(Air, fractionmass=1.);

//
// or use G4-NIST materials data base
//
G4NistManager* man = G4NistManager::Instance();
man->FindOrBuildMaterial("G4_SODIUM_IODIDE");

// print table
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the World
fDefaultMaterial  = vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeCalorParameters();
   
  //     
  // World
  //
  G4Box*
  solidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2); //its size
                         
  G4LogicalVolume*
  logicWorld = new G4LogicalVolume(solidWorld,              //its solid
                                   fDefaultMaterial,        //its material
                                   "World");                //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),         //at (0,0,0)
                                 logicWorld,                //its logical volume                                 
                                 "World",                   //its name
                                 0,                         //its mother  volume
                                 false,                     //no boolean operation
                                 0);                        //copy number
  
  //                               
  // Calorimeter
  //  
  G4Box* solidCalor=NULL; G4LogicalVolume* logicCalor=NULL; 
  G4Box* solidLayer=NULL; G4LogicalVolume* logicLayer=NULL;
  
  if (fCalorThickness > 0.)  
    { solidCalor = new G4Box("Calorimeter",                 //its name
                           fCalorThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);//size
                                 
      logicCalor = new G4LogicalVolume(solidCalor,             //its solid
                                             fDefaultMaterial, //its material
                                             "Calorimeter");   //its name
                                           
      new G4PVPlacement(0,                 //no rotation
                        G4ThreeVector(),   //at (0,0,0)
                        logicCalor,        //its logical volume
                        "Calorimeter",     //its name
                        logicWorld,        //its mother  volume
                        false,             //no boolean operation
                        0);                //copy number
  
  //                                 
  // Layer
  //
      solidLayer = new G4Box("Layer",                        //its name
                       fLayerThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); //size
                       
      logicLayer = new G4LogicalVolume(solidLayer,        //its solid
                                       fDefaultMaterial,  //its material
                                       "Layer");          //its name
      if (fNbOfLayers > 1)                                      
        new G4PVReplica("Layer",           //its name
                        logicLayer,        //its logical volume
                        logicCalor,        //its mother
                        kXAxis,            //axis of replication
                        fNbOfLayers,       //number of replica
                        fLayerThickness);  //width of replica
      else
        new G4PVPlacement(0,               //no rotation
                        G4ThreeVector(),   //at (0,0,0)
                        logicLayer,        //its logical volume                                 
                        "Layer",           //its name
                        logicCalor,        //its mother  volume
                        false,             //no boolean operation
                        0);                //copy number     
    }                                   
  
  //                               
  // Absorber
  //
  G4Box* solidAbsorber=NULL; G4LogicalVolume* logicAbsorber=NULL; 
  fPhysiAbsorber=NULL;  
  
  if (fAbsorberThickness > 0.) 
    { solidAbsorber = new G4Box("Absorber",                //its name
                          fAbsorberThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
      logicAbsorber = new G4LogicalVolume(solidAbsorber,        //its solid
                                          fAbsorberMaterial,    //its material
                                          fAbsorberMaterial->GetName()); //name
                                                
      fPhysiAbsorber = new G4PVPlacement(0,                  //no rotation
                           G4ThreeVector(-fGapThickness/2,0.,0.),  //its position
                                         logicAbsorber,     //its logical volume                    
                                         fAbsorberMaterial->GetName(), //its name
                                         logicLayer,        //its mother
                                         false,             //no boulean operat
                                         0);                //copy number
                                        
    }
  
  //                                 
  // Gap
  //
  G4Box* solidGap=NULL; G4LogicalVolume* logicGap=NULL; 
  fPhysiGap=NULL; 
  
  if (fGapThickness > 0.)
    { solidGap = new G4Box("Gap",
                           fGapThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
      logicGap = new G4LogicalVolume(solidGap,
                                     fGapMaterial,
                                     fGapMaterial->GetName());
                                           
      fPhysiGap = new G4PVPlacement(0,                    //no rotation
               G4ThreeVector(fAbsorberThickness/2,0.,0.), //its position
                             logicGap,                    //its logical volume               
                             fGapMaterial->GetName(),     //its name
                             logicLayer,                  //its mother
                             false,                       //no boulean operat
                             0);                          //copy number
    }
    
  PrintCalorParameters();     
  
  //                                        
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicCalor->SetVisAttributes(simpleBoxVisAtt);


  //
  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << fNbOfLayers << " layers of: [ "
         << fAbsorberThickness/mm << "mm of " << fAbsorberMaterial->GetName() 
         << " + "
         << fGapThickness/mm << "mm of " << fGapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) fAbsorberMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) fGapMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  fAbsorberThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapThickness(G4double val)
{
  // change Gap thickness and recompute the calorimeter parameters
  fGapThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fCalorSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int val)
{
  fNbOfLayers = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(fMagField) delete fMagField;           //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non nul
  { fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(fMagField);
    fieldMgr->CreateChordFinder(fMagField);
  } else {
    fMagField = 0;
    fieldMgr->SetDetectorField(fMagField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
