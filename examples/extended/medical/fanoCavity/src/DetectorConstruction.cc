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
// $Id: DetectorConstruction.cc,v 1.3 2007-10-08 12:05:02 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:pWall(0), pCavity(0)
{
  // default parameter values
  cavityThickness = 2*mm;
  cavityRadius    = 1*cm;      
  
  wallThickness = 5*mm;
  
  DefineMaterials();
  SetWallMaterial("Water");
  SetCavityMaterial("Water_vapor");
  
  // create commands for interactive definition of the detector  
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  //
  // define Elements
  //
  G4double z,a;
  
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);

  //
  // define materials
  //  
  G4Material* H2O = 
  new G4Material("Water", 1.0*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
  
  G4Material* vapor = 
  new G4Material("Water_vapor", 1.0*mg/cm3, 2);
  vapor->AddElement(H, 2);
  vapor->AddElement(O, 1);
  vapor->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
  
  G4Material* Air = 
  new G4Material("Air", 1.290*mg/cm3, 2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);
  
  new G4Material("Graphite",     6, 12.01*g/mole, 2.265*g/cm3);
  new G4Material("Graphite_gas", 6, 12.01*g/mole, 2.265*mg/cm3);  
  
  new G4Material("Aluminium",     13, 26.98*g/mole, 2.700*g/cm3);
  new G4Material("Aluminium_gas", 13, 26.98*g/mole, 2.700*mg/cm3);  
          
 G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
		   
  // Chamber
  //
  totalThickness = cavityThickness + 2*wallThickness;
  wallRadius     = cavityRadius + wallThickness;
  
  G4Tubs* 
  sChamber = new G4Tubs("Chamber",					//name
    		        0.,wallRadius,0.5*totalThickness,0.,twopi);	//size

  G4LogicalVolume*
  lChamber = new G4LogicalVolume(sChamber,		//solid
      			       wallMaterial,		//material
      			      "Chamber");		//name
				   
  pWall = new G4PVPlacement(0,				//no rotation
                             G4ThreeVector(),		//at (0,0,0)
                             lChamber,			//logical volume
                            "Wall",			//name
                             0,				//mother  volume
                             false,			//no boolean operation
                             0);			//copy number

  // Cavity
  //  			   
  G4Tubs*
  sCavity = new G4Tubs("Cavity",	
                       0.,cavityRadius,0.5*cavityThickness,0.,twopi);
		 
  G4LogicalVolume*		   			                      
  lCavity = new G4LogicalVolume(sCavity,		//shape
                                cavityMaterial,		//material
                                "Cavity");		//name
				
  pCavity = new G4PVPlacement(0,			//no rotation
                             G4ThreeVector(),		//at (0,0,0)
                             lCavity,			//logical volume
                            "Cavity",			//name
                             lChamber,			//mother  volume
                             false,			//no boolean operation
                             1);			//copy number
				
  PrintParameters();
    
  //
  //always return the root volume
  //  
  return pWall;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Wall is " << G4BestUnit(wallThickness,"Length")
         << " of " << wallMaterial->GetName() << " ( " 
	 << G4BestUnit(wallMaterial->GetDensity(),"Volumic Mass") << " )\n";
  G4cout << "     The Cavity is " << G4BestUnit(cavityThickness,"Length")
         << " of " << cavityMaterial->GetName() << " ( " 
	 << G4BestUnit(cavityMaterial->GetDensity(),"Volumic Mass") << " )";	 	 
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWallThickness(G4double value)
{
  wallThickness = value;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWallMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) wallMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityThickness(G4double value)
{
  cavityThickness = value;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityRadius(G4double value)
{
  cavityRadius  = value;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) cavityMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()
{
G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
