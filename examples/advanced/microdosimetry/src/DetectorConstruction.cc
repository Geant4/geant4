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
// -------------------------------------------------------------------
// $Id: DetectorConstruction.cc,v 1.2 2008/12/18 12:56:20 gunter Exp $
// -------------------------------------------------------------------

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
:physiWorld(NULL), logicWorld(NULL), solidWorld(NULL)
{}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{ 
  G4String name, symbol;             
  G4double density;            
  
  G4int ncomponents, natoms;
  
  // Elements  
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
 
  // Water 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="H2O"  , density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  
  // Default materials in setup.
  waterMaterial = H2O;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  // WORLD VOLUME
  
  WorldSizeX  = 1*mm; 
  WorldSizeY  = 1*mm; 
  WorldSizeZ  = 1*mm;

  solidWorld = new G4Box("World",				     //its name
			   WorldSizeX/2,WorldSizeY/2,WorldSizeZ/2);  //its size
  

  logicWorld = new G4LogicalVolume(solidWorld,	        //its solid
				   waterMaterial,	//its material
				   "World");		//its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

 
  // Visualization attributes
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  worldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(worldVisAtt);

  return physiWorld;
}
