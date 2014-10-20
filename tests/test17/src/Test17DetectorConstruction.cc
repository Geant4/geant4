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
//      ---------- Test17DetectorConstruction -------
//
//                by Vladimir Ivanchenko, 23 July 1999
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17DetectorConstruction.hh"
#include "Test17DetectorMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17DetectorConstruction::Test17DetectorConstruction():
  AbsorberMaterial(0),
  WorldMaterial(0),
  defaultWorld(true),
  solidWorld(0),
  logicWorld(0),
  physiWorld(0),
  solidAbsorber(0),
  logicAbsorber(0),
  physiAbsorber(0),
  magField(0)
{
  // default parameter values of the calorimeter
  AbsorberThickness = 100.0*mm;
  AbsorberSizeYZ    = 10.*cm;
  XposAbs           = 0.*cm ;
  NumberOfAbsorbers = 1;

  // create commands for interactive definition of the calorimeter
  detectorMessenger = new Test17DetectorMessenger(this);
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17DetectorConstruction::~Test17DetectorConstruction()
{
  delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::DefineMaterials()
{
  AbsorberMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  WorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Test17DetectorConstruction::Construct()
{
  // complete the Calor parameters definition and Print
  ComputeCalorParameters();

  //
  // World
  //
  solidWorld = new G4Box("World",				//its name
			 WorldSizeX/2.0,                        //its size X
                         WorldSizeYZ/2.0,WorldSizeYZ/2.0);      //its size YZ

  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   WorldMaterial,	//its material
                                   "World");		//its name

  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//no moving
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Absorber
  // 
  solidAbsorber = new G4Box("Absorber",	
                      AbsorberThickness/2,AbsorberSizeYZ/2,AbsorberSizeYZ/2); 
                          
  logicAbsorber = new G4LogicalVolume(solidAbsorber,    //its solid
    	                  AbsorberMaterial,             //its material
   	                 "Absorber");                   //its name
      			                  
  G4int copyNo=0;
  G4double x ;

  for (G4int j=0; j<NumberOfAbsorbers; j++)
  {
    x = XposAbs + AbsorberThickness * (G4double(j) + 0.5) ; 
    physiAbsorber = new G4PVPlacement(0,	   //no rotation
      		  G4ThreeVector(x,0.0,0.0),        //its position
                                "Absorber",        //its name
                                logicAbsorber,     //its logical volume
                                physiWorld,        //its mother
                                false,             //no boulean operat
                                copyNo);           //copy number
  }
  
  
  //
  //always return the physical World
  //
  PrintCalorParameters();  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The  WORLD   is made of "
         << G4BestUnit(WorldSizeX,"Length") << " of " 
	 << WorldMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the world is "
         << G4BestUnit(WorldSizeYZ,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " << NumberOfAbsorbers << " items of "
         <<G4BestUnit(AbsorberThickness,"Length")<< " of " 
	 << AbsorberMaterial->GetName();
  G4cout << ". The transverse size (YZ) is "
         << G4BestUnit(AbsorberSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the absorber " 
	 << G4BestUnit(XposAbs,"Length");
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetAbsorberMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if (pttoMaterial && pttoMaterial != AbsorberMaterial) {
    AbsorberMaterial = pttoMaterial;
    if(logicAbsorber) logicAbsorber->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if (pttoMaterial && pttoMaterial != WorldMaterial) {
    WorldMaterial = pttoMaterial;
    if(logicWorld) logicWorld->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetNumberOfAbsorbers(G4int val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  NumberOfAbsorbers = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetAbsorberSizeYZ(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  AbsorberSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetWorldSizeX(G4double val)
{
  WorldSizeX = val;
  defaultWorld = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  WorldSizeYZ = val;
  defaultWorld = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetAbsorberXpos(G4double val)
{
  XposAbs  = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;		//delete the existing magn field

  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
