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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
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

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  for(G4int i=0; i<kMaxAbsor; ++i) { 
    fAbsorMaterial[i] = nullptr; 
    fAbsorThickness[i] = 0.0;
    fSolidAbsor[i] = nullptr;
    fLogicAbsor[i] = nullptr;
    fPhysiAbsor[i] = nullptr;
  } 

  // default parameter values of the calorimeter
  fNbOfAbsor = 2;
  fAbsorThickness[1] = 36*mm;
  fAbsorThickness[2] = 4*mm;
  fNbOfLayers        = 50;
  fCalorSizeYZ       = 1.5*m;
  ComputeCalorParameters();

  // materials
  DefineMaterials();
  SetWorldMaterial("Galactic");
  SetAbsorMaterial(1,"Iron");
  SetAbsorMaterial(2,"Scintillator");

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // This function illustrates the possible ways to define materials using 
  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);
  //
  // define Elements
  //
  G4Element* H  = manager->FindOrBuildElement(1);
  G4Element* C  = manager->FindOrBuildElement(6);
  G4Element* O  = manager->FindOrBuildElement(8);
  //
  // define an Element from isotopes, by relative abundance
  //
  G4int iz, n;                       //iz=number of protons  in an isotope;
                                     // n=number of nucleons in an isotope;                            
  G4int   ncomponents;
  G4double z, a;                                     
  G4double abundance;                                     

  G4Isotope* U5 = new G4Isotope("U235", iz=92, n=235, a=235.01*g/mole);
  G4Isotope* U8 = new G4Isotope("U238", iz=92, n=238, a=238.03*g/mole);

  G4Element* U  = new G4Element("enriched Uranium", "U", ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);

  //
  // define simple materials
  //
  G4double density;

  new G4Material("liquidH2",    z=1.,  a= 1.008*g/mole,  density= 70.8*mg/cm3);
  new G4Material("Aluminium",   z=13., a= 26.98*g/mole,  density= 2.700*g/cm3);
  new G4Material("liquidArgon", z=18,  a= 39.948*g/mole, density= 1.396*g/cm3);
  new G4Material("Titanium",    z=22., a= 47.867*g/mole, density= 4.54*g/cm3);
  new G4Material("Iron",        z=26., a= 55.85*g/mole,  density= 7.870*g/cm3);
  new G4Material("Copper",      z=29., a= 63.55*g/mole,  density= 8.960*g/cm3);
  new G4Material("Tungsten",    z=74., a= 183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",        z=79., a= 196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",        z=82., a= 207.20*g/mole, density= 11.35*g/cm3);  
  new G4Material("Uranium",     z=92., a= 238.03*g/mole, density= 18.95*g/cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //
  G4int natoms;

  G4Material* H2O = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  H2O->SetChemicalFormula("H_2O");
  
  G4Material* CH = 
  new G4Material("Polystyrene", density= 1.032*g/cm3, ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* Sci = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);
  
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  //
  // examples of gas in non STP conditions
  //
  G4double temperature, pressure;
  
  G4Material* CO2 = 
  new G4Material("CarbonicGas", density= 27.*mg/cm3, ncomponents=2,
                 kStateGas, temperature= 325.*kelvin, pressure= 50.*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);
  
  new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);
  //
  // example of vacuum
  //
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                             kStateGas,temperature,pressure);

  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fLayerThickness = 0.;
  for (G4int iAbs=1; iAbs<=fNbOfAbsor; iAbs++) {
    fLayerThickness += fAbsorThickness[iAbs];
  }
  fCalorThickness = fNbOfLayers*fLayerThickness;     
  fWorldSizeX = 1.2*fCalorThickness; 
  fWorldSizeYZ = 1.2*fCalorSizeYZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if(fPhysiWorld) { return fPhysiWorld; }
  // complete the Calor parameters definition
  ComputeCalorParameters();

  //
  // World
  //
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);  //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,        //its solid
                                    fWorldMaterial,     //its material
                                    "World");           //its name

  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                                  G4ThreeVector(),      //at (0,0,0)
                                  fLogicWorld,          //its fLogical volume
                                  "World",              //its name
                                  0,                    //its mother  volume
                                  false,                //no boolean operation
                                  0);                   //copy number
  //
  // Calorimeter
  //

  fSolidCalor = new G4Box("Calorimeter",                                  
                       fCalorThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);

  fLogicCalor = new G4LogicalVolume(fSolidCalor,               
                                    fWorldMaterial,     
                                    "Calorimeter");      

  fPhysiCalor = new G4PVPlacement(0,                     //no rotation
                                 G4ThreeVector(),        //at (0,0,0)
                                 fLogicCalor,            //its fLogical volume
                                 "Calorimeter",          //its name
                                 fLogicWorld,            //its mother  volume
                                 false,                  //no boolean operation
                                 0);                     //copy number

  //
  // Layers
  //

  fSolidLayer = new G4Box("Layer",                               
                          fLayerThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);

  fLogicLayer = new G4LogicalVolume(fSolidLayer,      
                                    fWorldMaterial,   
                                    "Layer");              
  if (fNbOfLayers > 1) {
    fPhysiLayer = new G4PVReplica("Layer",              
                                  fLogicLayer,     
                                  fLogicCalor,      
                                  kXAxis,              
                                  fNbOfLayers,            
                                  fLayerThickness);     
  } else {
    fPhysiLayer = new G4PVPlacement(0,                   
                                   G4ThreeVector(),     
                                   fLogicLayer,           
                                   "Layer",             
                                   fLogicCalor,         
                                   false,             
                                   0);                    
  }
  //
  // Absorbers
  //

  G4double xfront = -0.5*fLayerThickness;
  for (G4int k=1; k<=fNbOfAbsor; ++k) {
    fSolidAbsor[k] = new G4Box("Absorber",                //its name
                         fAbsorThickness[k]/2,fCalorSizeYZ/2,fCalorSizeYZ/2);

    fLogicAbsor[k] = new G4LogicalVolume(fSolidAbsor[k],    //its solid
                                         fAbsorMaterial[k], //its material
                                         fAbsorMaterial[k]->GetName());

    G4double xcenter = xfront+0.5*fAbsorThickness[k];
    xfront += fAbsorThickness[k];
    fPhysiAbsor[k] = new G4PVPlacement(0,              
                         G4ThreeVector(xcenter,0.,0.),
                         fLogicAbsor[k],               
                         fAbsorMaterial[k]->GetName(),
                         fLogicLayer,                  
                         false,                      
                         k);                                //copy number

  }

  PrintCalorParameters();

  //always return the fPhysical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4int prec = 4, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  G4double totLength(0.), totRadl(0.), totNuclear(0.);
  
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << fNbOfLayers << " layers of:";
  for (G4int i=1; i<=fNbOfAbsor; ++i) {
    G4Material* material = fAbsorMaterial[i];
    G4double radl = material->GetRadlen();
    G4double nuclearl = material->GetNuclearInterLength();
    G4double sumThickness = fNbOfLayers*fAbsorThickness[i];
    G4double nbRadl = sumThickness/radl;
    G4double nbNuclearl = sumThickness/nuclearl;
    totLength += sumThickness;
    totRadl += nbRadl;
    totNuclear += nbNuclearl;
    G4cout << "\n   " << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
           << std::setw(wid) << G4BestUnit(fAbsorThickness[i],"Length")
           << "  --->  sum = " << std::setw(wid) << G4BestUnit(sumThickness,"Length")
           << " = " << std::setw(wid) << nbRadl << " Radl "
           << " = " << std::setw(wid) << nbNuclearl << " NuclearInteractionLength " ;
  }
  G4cout << "\n\n                       total thickness = " 
                  << std::setw(wid) << G4BestUnit(totLength,"Length")
         << " = " << std::setw(wid)<< totRadl << " Radl "
         << " = " << std::setw(wid)<< totNuclear << " NuclearInteractionLength "
         << G4endl;
         
  G4cout << "                     transverse sizeYZ = " 
                  << std::setw(wid) << G4BestUnit(fCalorSizeYZ,"Length")
         << G4endl;         
  G4cout << "-------------------------------------------------------------\n";
  
  G4cout << "\n" << fWorldMaterial << G4endl;    
  for (G4int j=1; j<=fNbOfAbsor; ++j) {
    G4cout << "\n" << fAbsorMaterial[j] << G4endl;
  }
  G4cout << "\n-------------------------------------------------------------\n";
  
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& material)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if(pttoMaterial) { 
    fWorldMaterial = pttoMaterial;
    if(fLogicWorld) {
      fLogicWorld->SetMaterial(fWorldMaterial);
      fLogicLayer->SetMaterial(fWorldMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int ival)
{
  // set the number of Layers
  //
  if (ival < 1)
    { G4cout << "\n --->warning from SetfNbOfLayers: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fNbOfLayers = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << kMaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  fNbOfAbsor = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int ival,
                                            const G4String& material)
{
  // search the material by its name
  //
  if (ival > fNbOfAbsor || ival <= 0)
    { G4cout << "\n --->warning from SetAbsorMaterial: absor number "
             << ival << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
    fAbsorMaterial[ival] = pttoMaterial;
    if(fLogicAbsor[ival]) {
      fLogicAbsor[ival]->SetMaterial(pttoMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();    
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int ival, G4double val)
{
  // change Absorber thickness
  //
  if (ival > fNbOfAbsor || ival <= 0)
    { G4cout << "\n --->warning from SetAbsorThickness: absor number "
             << ival << " out of range. Command refused" << G4endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorThickness[ival] = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfCalorSizeYZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fCalorSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

void DetectorConstruction::ConstructSDandField()
{
  if ( fFieldMessenger.Get() == nullptr ) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    G4GlobalMagFieldMessenger* msg =
      new G4GlobalMagFieldMessenger(fieldValue);
    //msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put( msg );        
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
