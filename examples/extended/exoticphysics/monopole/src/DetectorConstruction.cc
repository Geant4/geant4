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
// $Id: DetectorConstruction.cc,v 1.2 2010-06-04 19:03:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MonopoleFieldSetup.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  // default parameter values
  absorSizeX = absorSizeYZ = 10 * cm;
  worldSizeX = worldSizeYZ = 1.2 * absorSizeX;
  maxStepSize = 5 * mm;

  worldMaterial = absorMaterial = 0;
  magField = 0;
  lAbsor   = 0;
  fMFieldSetup = 0;

  DefineMaterials();
  SetMaterial("G4_Al");

  // create commands for interactive definition of the detector
  detectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  return ConstructVolumes();
}

void DetectorConstruction::DefineMaterials()
{ 
  /***********************     define Elements     *********************/
  G4double z, a;

  G4Element* N = new G4Element("Nitrogen", "N", z= 7, a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z= 8, a= 16.00*g/mole);

  /***********************    define Materials    **********************/
  G4double density, temperature, pressure;
  G4int    ncomponents;
  G4double fractionmass;
 
  G4Material* Air = new G4Material("Air"  , density = 1.290 * mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18 * pascal;
  temperature = 2.73 * kelvin;

  G4Material* vacuum = new G4Material("Galactic", z = 1, a = 1.008 * g/mole, density,
                          kStateGas, temperature, pressure);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //default materials
  worldMaterial = vacuum;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  /****************************    World   *****************************/
  G4Box * sWorld = new G4Box("world",					//name
	       worldSizeX / 2, worldSizeYZ / 2, worldSizeYZ / 2);	//dimensions

  G4LogicalVolume * lWorld = new G4LogicalVolume(sWorld,		//shape
						worldMaterial,		//material
						"world");		//name

  G4VPhysicalVolume * pWorld = new G4PVPlacement(0,	//no rotation
  			   G4ThreeVector(),		//at (0,0,0)
                           lWorld,			//logical volume
			   "world",			//name
                           0,	       		        //mother  volume
                           false,			//no boolean operation
                           0);				//copy number


  /**************************    Absorber    ***************************/
  G4Box * sAbsor = new G4Box("Absorber",			//name
                   absorSizeX / 2, absorSizeYZ / 2, absorSizeYZ / 2);	//dimensions
		   			                      
  lAbsor = new G4LogicalVolume(sAbsor,			//shape
                               absorMaterial,		//material
                              "Absorber");		//name
  
                              
           new G4PVPlacement(0,				//no rotation
  			   G4ThreeVector(),		//at (0,0,0)
                           lAbsor,			//logical volume
			   "Absorber",			//name
                           lWorld,	       		//mother  volume
                           false,			//no boolean operation
                           0);				//copy number
  lAbsor->SetUserLimits(new G4UserLimits(maxStepSize));

  PrintParameters();

  /************     always return the World volume     *****************/
  return pWorld;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(absorSizeX, "Length")
         << " of " << absorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSizeX(G4double value)
{
  absorSizeX = value; worldSizeX = 1.2 * absorSizeX;
}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSizeYZ(G4double value)
{
  absorSizeYZ = value; 
  worldSizeYZ = 1.2 * absorSizeYZ;
}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material * pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial) absorMaterial = pttoMaterial;
  else  G4cout << "Material does not exist in DB" << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager * fieldMgr = 
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if (magField) { delete magField; }	//delete the existing magn field

  //fMFieldSetup = G4MonopoleFieldSetup::GetMonopoleFieldSetup(); // create the field

  
  if (fieldValue != 0.)			// create a new one if non nul
    {
      magField = new G4UniformMagField(G4ThreeVector(0., 0., fieldValue));        
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
    }
   else
    {
      magField = 0;
      fieldMgr->SetDetectorField(magField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStepSize(G4double step_)
{
  maxStepSize = step_;
  if(lAbsor) lAbsor->SetUserLimits(new G4UserLimits(maxStepSize));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManager.hh" 

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
