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
// $Id: Tst20DetectorConstruction.cc,v 1.9 2007-11-09 18:33:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#include "Tst20DetectorConstruction.hh"
#include "Tst20DetectorMessenger.hh"

#include "Tst20CalorimeterSD.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4ios.hh"


Tst20DetectorConstruction::Tst20DetectorConstruction(): worldChanged(false),
							absorberMaterial(0),
							worldMaterial(0),
							solidWorld(0),
							logicWorld(0),
							physiWorld(0),
							solidAbsorber(0),
							logicAbsorber(0),
							physiAbsorber(0),
							magneticField(0),
							calorimeterSD(0)
{
  // Default parameter values of the calorimeter
  worldSizeZ = 5*micrometer;
  worldSizeR = 5*micrometer;
  absorberThickness = 1*micrometer;
  absorberRadius   = 1*micrometer;
  zAbsorber = 0.*cm ;

  // Create commands for interactive definition of the calorimeter  
  detectorMessenger = new Tst20DetectorMessenger(this);

  // Define all necessary materials (not deleted!)
  DefineMaterials();
}


Tst20DetectorConstruction::~Tst20DetectorConstruction()
{ 
  delete detectorMessenger;
}


G4VPhysicalVolume* Tst20DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}


void Tst20DetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials
 
  G4String name;
  G4String symbol;   
  G4double a; // mass of a mole
  G4double z; // mean number of protons

  // define Elements

  a = 1.01*g/mole;
  z = 1.;
  symbol="H";
  name="Hydrogen";
  G4Element* elH  = new G4Element(name, symbol, z, a);

  a = 14.01*g/mole;
  name = "Nitrogen"; 
  symbol = "N";
  z = 7.;
  G4Element* elN  = new G4Element(name, symbol, z, a);

  a = 16.00*g/mole;
  name = "Oxygen";
  symbol = "O";
  z = 8.;
  G4Element* elO  = new G4Element(name, symbol, z, a);

  G4double density;             
  G4int nComponents;
  G4int nAtoms;
  G4double fractionMass;

  // Define a material from elements.   case 1: chemical molecule

  density = 1.000*g/cm3;
  name = "Water";
  nComponents = 2;
  nAtoms = 2;
  G4Material* water = new G4Material(name, density, nComponents);
  nAtoms = 2;
  water->AddElement(elH, nAtoms);
  nAtoms = 1;
  water->AddElement(elO, nAtoms);

  // Define a material from elements.   case 2: mixture by fractional mass

  density = 1.290*mg/cm3;
  name = "Air";
  nComponents = 2;
  G4Material* air = new G4Material(name, density, nComponents);
  fractionMass = 0.7;
  air->AddElement(elN, fractionMass);
  fractionMass = 0.3;
  air->AddElement(elO, fractionMass);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // Default materials of the calorimeter and world
  absorberMaterial = water;
  worldMaterial = air;
}


G4VPhysicalVolume* Tst20DetectorConstruction::ConstructCalorimeter()
{
  ComputeCalorParameters();
  PrintCalorParameters();
  CleanGeometry();
     
  // World

  solidWorld = new G4Tubs("World",				//its name
			  0.,worldSizeR,worldSizeZ/2.,0.,twopi)       ;//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   absorberMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(0,0,0),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
                           
  // Absorber
 
  if (absorberThickness > 0.) 
    {
      solidAbsorber = new G4Tubs("Absorber",0.,absorberRadius,absorberThickness/2.,0.,twopi); 
                          
      logicAbsorber = new G4LogicalVolume(solidAbsorber,    //its solid
      			                  absorberMaterial, //its material
      			                  "Absorber");      //its name
      			                  
      physiAbsorber = new G4PVPlacement(0,		   //no rotation
					G4ThreeVector(0.,0.,zAbsorber),  //its position
                                        "Absorber",        //its name
                                        logicAbsorber,     //its logical volume
                                        physiWorld,        //its mother
                                        false,             //no boulean operat
                                        0);                //copy number
                                        
    }
                               
  // Sensitive Detectors: Absorber 
 
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();

  if (!calorimeterSD)
    {
      calorimeterSD = new Tst20CalorimeterSD("CalorSD",this);
      sdManager->AddNewDetector( calorimeterSD );
    }
  if (logicAbsorber) logicAbsorber->SetSensitiveDetector(calorimeterSD);

  return physiWorld;
}


void Tst20DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The  WORLD   is made of " 
	 << worldSizeZ/mm << "mm of " << worldMaterial->GetName()
	 << ", the transverse size (R) of the world is " << worldSizeR/mm << " mm " << G4endl;
  G4cout << " The ABSORBER is made of " 
	 << absorberThickness/mm << "mm of " << absorberMaterial->GetName() 
	 << ", the transverse size (R) is " << absorberRadius/mm << " mm " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " << zAbsorber/mm << "  mm" << G4endl;
}


void Tst20DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  // Search the material by its name   
  G4Material* material;
  for (size_t j=0 ; j<G4Material::GetNumberOfMaterials() ; j++)
    { 
      material = (*materialTable)[j];     
      if (material->GetName() == materialChoice)
	{
	  absorberMaterial = material;
	  logicAbsorber->SetMaterial(material); 
	  return;
	}             
    }
  G4cout << "Unvalid material" << G4endl;
}


void Tst20DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // Get the pointer to the material table
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  // Search the material by its name   
  G4Material* material;
  for (size_t j=0 ; j<G4Material::GetNumberOfMaterials() ; j++)
    { 
      material = (*materialTable)[j];     
      if (material->GetName() == materialChoice)
	{
	  worldMaterial = material;
	  logicWorld->SetMaterial(material); 
	  return;
	}             
    }
  G4cout << "Unvalid material" << G4endl;
}


void Tst20DetectorConstruction::SetAbsorberThickness(G4double value)
{
  // Change Absorber thickness and recompute the calorimeter parameters
  absorberThickness = value;
  ComputeCalorParameters();
}  


void Tst20DetectorConstruction::SetAbsorberRadius(G4double value)
{
  // Change the transverse size and recompute the calorimeter parameters
  absorberRadius = value;
  ComputeCalorParameters();
}  


void Tst20DetectorConstruction::SetWorldSizeZ(G4double value)
{
  worldChanged = true;
  worldSizeZ = value;
  ComputeCalorParameters();
}  


void Tst20DetectorConstruction::SetWorldSizeR(G4double value)
{
  worldChanged=true;
  worldSizeR = value;
  ComputeCalorParameters();
}  


void Tst20DetectorConstruction::SetAbsorberZpos(G4double value)
{
  zAbsorber  = value;
  ComputeCalorParameters();
}  


void Tst20DetectorConstruction::SetMagField(G4double fieldValue)
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if (magneticField) delete magneticField;		//delete the existing magnetic field
  
  if (fieldValue != 0.)			// create a new one if non null
    { 
      magneticField = new G4UniformMagField(G4ThreeVector(fieldValue,0.,0.));        
      fieldManager->SetDetectorField(magneticField);
      fieldManager->CreateChordFinder(magneticField);
    } 
  else 
    {
      magneticField = 0;
      fieldManager->SetDetectorField(magneticField);
    }
}

  
void Tst20DetectorConstruction::CleanGeometry()
{
  if (physiWorld)
    {
      G4PhysicalVolumeStore::Clean();
      G4LogicalVolumeStore::Clean();
      G4SolidStore::Clean();
    }
}

  
void Tst20DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}


void Tst20DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     if (! worldChanged)
     {
       worldSizeR = 2. * absorberRadius ;
       worldSizeZ = 2. * absorberThickness ;
     }
     
     zStartAbs = zAbsorber - 0.5 * absorberThickness; 
     zEndAbs = zAbsorber + 0.5 * absorberThickness; 

}
