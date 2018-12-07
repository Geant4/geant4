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
/// \file electromagnetic/TestEm3/src/DetectorConstruction.cc
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
 :G4VUserDetectorConstruction(),
  fWorldMaterial(nullptr),fSolidWorld(nullptr),fLogicWorld(nullptr),
  fPhysiWorld(nullptr),fSolidCalor(nullptr),fLogicCalor(nullptr),
  fPhysiCalor(nullptr),fSolidLayer(nullptr),fLogicLayer(nullptr),
  fPhysiLayer(nullptr)
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
  fAbsorThickness[1] = 2.3*mm;
  fAbsorThickness[2] = 5.7*mm;
  fNbOfLayers        = 50;
  fCalorSizeYZ       = 40.*cm;
  ComputeCalorParameters();

  // materials
  DefineMaterials();
  SetWorldMaterial("Galactic");
  SetAbsorMaterial(1,"G4_Pb");
  SetAbsorMaterial(2,"G4_lAr");

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
  G4double z,a;

  G4Element* H  = manager->FindOrBuildElement(1);
  G4Element* C  = manager->FindOrBuildElement(6);
  G4Element* N  = manager->FindOrBuildElement(7);
  G4Element* O  = manager->FindOrBuildElement(8);
  G4Element* Si = manager->FindOrBuildElement(14);
  G4Element* Ge = manager->FindOrBuildElement(32);
  G4Element* Sb = manager->FindOrBuildElement(51);
  G4Element* I  = manager->FindOrBuildElement(53);
  G4Element* Cs = manager->FindOrBuildElement(55);
  G4Element* Pb = manager->FindOrBuildElement(82);
  G4Element* Bi = manager->FindOrBuildElement(83);

  //
  // define an Element from isotopes, by relative abundance
  //
  G4int iz, n;                       //iz=number of protons  in an isotope;
                                     // n=number of nucleons in an isotope;
  G4int   ncomponents;                                     
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
  new G4Material("Titanium",    z=22., a= 47.867*g/mole, density= 4.54*g/cm3);
  new G4Material("Iron",        z=26., a= 55.85*g/mole,  density= 7.870*g/cm3);
  new G4Material("Copper",      z=29., a= 63.55*g/mole,  density= 8.960*g/cm3);
  new G4Material("Tungsten",    z=74., a= 183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",        z=79., a= 196.97*g/mole, density= 19.32*g/cm3);
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

  G4Material* Lct =
  new G4Material("Lucite", density= 1.185*g/cm3, ncomponents=3);
  Lct->AddElement(C, 59.97*perCent);
  Lct->AddElement(H, 8.07*perCent);
  Lct->AddElement(O, 31.96*perCent);

  G4Material* Sili = 
  new G4Material("Silicon", density= 2.330*g/cm3, ncomponents=1);
  Sili->AddElement(Si, natoms=1);

  G4Material* SiO2 = 
  new G4Material("quartz", density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  G4Material* G10 = 
  new G4Material("NemaG10", density= 1.700*g/cm3, ncomponents=4);
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);

  G4Material* CsI = 
  new G4Material("CsI", density= 4.534*g/cm3, ncomponents=2);
  CsI->AddElement(Cs, natoms=1);
  CsI->AddElement(I , natoms=1);
  CsI->GetIonisation()->SetMeanExcitationEnergy(553.1*eV);

  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);

  //SiNx
  density= 3.1 *g/cm3;
  G4Material* SiNx= new G4Material("SiNx", density, ncomponents=3);
  SiNx-> AddElement(Si, 300);
  SiNx-> AddElement(N, 310);
  SiNx-> AddElement(H, 6);

  //
  // define gaseous materials using G4 NIST database 
  //
  G4double fractionmass;
  
  G4Material* Air = manager->FindOrBuildMaterial("G4_AIR");
  manager->ConstructNewGasMaterial("Air20","G4_AIR",293.*kelvin,1.*atmosphere);

  G4Material* lAr = manager->FindOrBuildMaterial("G4_lAr");
  G4Material* lArEm3 = new G4Material("liquidArgon", density= 1.390*g/cm3,
                                                                ncomponents=1);
  lArEm3->AddMaterial(lAr, fractionmass=1.0);

  //
  // define a material from elements and others materials (mixture of mixtures)
  //

  G4Material* Lead = new G4Material("Lead",density=11.35*g/cm3,ncomponents=1);
  Lead->AddElement(Pb, fractionmass=1.0);

  G4Material* LeadSb = new G4Material("LeadSb", density=11.35*g/cm3, 
                                                                ncomponents=2);
  LeadSb->AddElement(Sb, fractionmass=4.*perCent);
  LeadSb->AddElement(Pb, fractionmass=96.*perCent);

  G4Material* Aerog = new G4Material("Aerogel", density= 0.200*g/cm3,
                                                                ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
  Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
  Aerog->AddElement (C   , fractionmass= 0.1*perCent);

  //
  // examples of gas in non STP conditions
  //
  G4double temperature, pressure;
  
  G4Material* CO2 = 
  new G4Material("CarbonicGas", density= 27.*mg/cm3, ncomponents=2,
                 kStateGas, temperature= 325.*kelvin, pressure= 50.*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* steam = 
  new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1,
                  kStateGas, temperature= 273*kelvin, pressure= 1*atmosphere);
  steam->AddMaterial(H2O, fractionmass=1.);
  
  new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);
  //
  // examples of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                             kStateGas,temperature,pressure);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;         //from PhysicalConstants.h
  G4Material* beam = 
  new G4Material("Beam", density, ncomponents=1,
                         kStateGas,temperature,pressure);
  beam->AddMaterial(Air, fractionmass=1.);

  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
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
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << fNbOfLayers << " layers of:";
  for (G4int i=1; i<=fNbOfAbsor; ++i) {
    G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
           << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length");
  }
  G4cout << "\n-------------------------------------------------------------\n";
  
  G4cout << "\n" << fWorldMaterial << G4endl;    
  for (G4int j=1; j<=fNbOfAbsor; ++j) {
    G4cout << "\n" << fAbsorMaterial[j] << G4endl;
  }
  G4cout << "\n-------------------------------------------------------------\n";
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
