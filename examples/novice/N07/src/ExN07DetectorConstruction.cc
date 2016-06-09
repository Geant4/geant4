//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN07DetectorConstruction.cc,v 1.3 2003/04/08 15:47:00 asaim Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

#include "ExN07DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "ExN07DetectorMessenger.hh"
#include "ExN07CalorimeterSD.hh"
#include "ExN07GapParameterisation.hh"
#include "ExN07PrimaryGeneratorAction.hh"

ExN07DetectorConstruction::ExN07DetectorConstruction()
:numberOfLayers(40),worldMaterial(0),absorberMaterial(0),gapMaterial(0),
 gapSolid(0),worldLogical(0),worldPhysical(0),serial(false)
{
  for(size_t i=0;i<3;i++)
  {
    calorLogical[i] = 0;
    layerLogical[i] = 0;
    calorPhysical[i] = 0;
    layerPhysical[i] = 0;
  }

  gapParam = new ExN07GapParameterisation;
  gapParam->SetNumberOfLayers(numberOfLayers);
  ExN07CalorimeterSD::SetNumberOfLayers(numberOfLayers);
  DefineMaterials();
  detectorMessenger = new ExN07DetectorMessenger(this);
}

ExN07DetectorConstruction::~ExN07DetectorConstruction()
{ delete detectorMessenger;}


void ExN07DetectorConstruction::DefineMaterials()
{ 
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  
  G4int iz, n;                       //iz=number of protons  in an isotope; 
                                     // n=number of nucleons in an isotope;

  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //
  // define an Element from isotopes, by relative abundance 
  //

  G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
  G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

  G4Element* U  = new G4Element(name="enriched Uranium",symbol="U",ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);

  //
  // define simple materials
  //

  new G4Material(name="Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);

  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

  //
  // define a material from elements.   case 1: chemical molecule
  //
 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  //
  // examples of vacuum
  //

  density     = universe_mean_density;
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole,
                                    density,kStateGas,temperature,pressure);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //default materials of the calorimeter
  worldMaterial    = Vacuum;
  absorberMaterial = Pb;
  gapMaterial      = lAr;
  gapParam->SetAbsorberMaterial(absorberMaterial);
  gapParam->SetGapMaterial(gapMaterial);
}

G4VPhysicalVolume* ExN07DetectorConstruction::Construct()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",2.*m,2.*m,4.*m);
  worldLogical = new G4LogicalVolume(worldSolid,worldMaterial,"World");
  worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",
                        0,false,0);
  
  //                               
  // Calorimeter
  //  
  G4VSolid* calorSolid = new G4Box("Calor",0.5*m,0.5*m,1.*m);
  G4int i;
  G4String calNam[3] = {"Cal-A","Cal-B","Cal-C"};
  for(i=0;i<3;i++)
  {
    calorLogical[i] = new G4LogicalVolume(calorSolid,absorberMaterial,calNam[i]);
    if(serial)
    {
      calorPhysical[i] = new G4PVPlacement(0,
                 G4ThreeVector(0.,0.,G4double(i-1)*2.*m),
                 calorLogical[i],calNam[i],worldLogical,false,i);
    }
    else
    {
      calorPhysical[i] = new G4PVPlacement(0,
                 G4ThreeVector(0.,G4double(i-1)*m,0.),
                 calorLogical[i],calNam[i],worldLogical,false,i);
    }
  }
 
  //                                 
  // Layers -- material is parameterised
  //
  gapSolid = new G4Box("Gap",0.5*m,0.5*m,1.*m/G4double(2*numberOfLayers));
  for(i=0;i<3;i++)
  {
    layerLogical[i] = new G4LogicalVolume(gapSolid,absorberMaterial,"Layer");
    layerPhysical[i] = 
       new G4PVParameterised("Layer",layerLogical[i],calorLogical[i],kZAxis,
                          2*numberOfLayers,gapParam);
  }
   
  //
  // Regions
  //
  G4String regName[] = {"Calor-A","Calor-B","Calor-C"};
  for(i=0;i<3;i++)
  {
    G4Region* aRegion = new G4Region(regName[i]);
    calorLogical[i]->SetRegion(aRegion);
    aRegion->AddRootLogicalVolume(calorLogical[i]);
  }

  //                               
  // Sensitive Detectors: Absorber and Gap
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String detName[] = {"CalorSD-A","CalorSD-B","CalorSD-C"};
  for(i=0;i<3;i++)
  {
    G4VSensitiveDetector* calorSD = new ExN07CalorimeterSD(detName[i]);
    SDman->AddNewDetector(calorSD);
    layerLogical[i]->SetSensitiveDetector(calorSD);
  }
  
  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  for(i=0;i<3;i++)
  { 
    calorLogical[i]->SetVisAttributes(simpleBoxVisAtt);
    layerLogical[i]->SetVisAttributes(simpleBoxVisAtt);
  }
  
  PrintCalorParameters();
  return worldPhysical;
}

void ExN07DetectorConstruction::PrintCalorParameters() const
{
  G4cout << "--------------------------------------------------------" << G4endl;
  if(serial)
  { G4cout << " Calorimeters are placed in serial." << G4endl; }
  else
  { G4cout << " Calorimeters are placed in parallel." << G4endl; }
  G4cout << " Absorber is made of " << absorberMaterial->GetName() << G4endl;
  G4cout << " Gap is made of " << gapMaterial->GetName() << G4endl;
  G4cout << "--------------------------------------------------------" << G4endl;
}

void ExN07DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if(pttoMaterial)
  {
    absorberMaterial = pttoMaterial;
    gapParam->SetAbsorberMaterial(pttoMaterial);
    for(size_t i=0;i<3;i++)
    {
      calorLogical[i]->SetMaterial(absorberMaterial);
      layerLogical[i]->SetMaterial(absorberMaterial);
    }
  }
  else
  { G4cerr << materialChoice << " is not defined. - Command is ignored." << G4endl; }
}

G4String ExN07DetectorConstruction::GetAbsorberMaterial() const
{ return absorberMaterial->GetName(); }

void ExN07DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if(pttoMaterial)
  {
    gapMaterial = pttoMaterial;
    gapParam->SetGapMaterial(pttoMaterial);
  }
  else
  { G4cerr << materialChoice << " is not defined. - Command is ignored." << G4endl; }
}

G4String ExN07DetectorConstruction::GetGapMaterial() const
{ return gapMaterial->GetName(); }

void ExN07DetectorConstruction::SetSerialGeometry(G4bool ser)
{
  if(serial==ser) return;
  serial=ser;
  for(G4int i=0;i<3;i++)
  {
    if(serial)
    { calorPhysical[i]->SetTranslation(G4ThreeVector(0.,0.,G4double(i-1)*2.*m)); }
    else
    { calorPhysical[i]->SetTranslation(G4ThreeVector(0.,G4double(i-1)*m,0.)); }
  }
  ExN07PrimaryGeneratorAction* gen = (ExN07PrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  gen->SetSerial(serial);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ExN07DetectorConstruction::SetNumberOfLayers(G4int nl)
{
  numberOfLayers = nl;
  gapSolid->SetZHalfLength(1.*m/G4double(2*numberOfLayers));
  gapParam->SetNumberOfLayers(nl);
  for(size_t i=0;i<3;i++)
  { 
    if(layerPhysical[i])
    {
      calorLogical[i]->RemoveDaughter(layerPhysical[i]);
      delete layerPhysical[i];
    }
    layerPhysical[i] = 
       new G4PVParameterised("Layer",layerLogical[i],calorLogical[i],kZAxis,
                          2*numberOfLayers,gapParam);
  }
  ExN07CalorimeterSD::SetNumberOfLayers(nl);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

