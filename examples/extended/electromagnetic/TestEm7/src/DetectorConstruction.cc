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
// $Id: DetectorConstruction.cc,v 1.10 2008-04-21 13:13:30 vnivanch Exp $
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
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  // default parameter values
  absorSizeX = absorSizeYZ = 20*cm;
  worldSizeX = worldSizeYZ = 1.2*absorSizeX;
  
  worldMaterial = absorMaterial = 0;
  magField = 0;
  lAbsor   = 0;
  
  tallySize     = G4ThreeVector();
  tallyMaterial = 0;
  tallyMass     = 0.; 
  tallyNumber   = 0;
  tallyPosition = new G4ThreeVector[MaxTally];
  lTally        = 0;
  
  DefineMaterials();
  SetMaterial("Water");
  SetTallyMaterial("Water");

  // create commands for interactive definition of the detector  
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete [] tallyPosition; delete detectorMessenger;}

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
  G4double z, a;

  G4Element* H = new G4Element("Hydrogen", "H", z= 1, a= 1.008*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z= 7, a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z= 8, a= 16.00*g/mole);

  //
  // define Materials.
  //
  G4double density, temperature, pressure;
  G4int    ncomponents, natoms;
  G4double fractionmass;
 
  G4Material* H2O = 
    new G4Material("Water", density= 1.0*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  G4Material* Air = 
    new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* vac = new G4Material( "TechVacuum", density, 1,
                           kStateGas, temperature, pressure );
  vac->AddMaterial( Air, 1. );

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* vacuum = 
    new G4Material("Galactic",z= 1,a= 1.008*g/mole,density,
		   kStateGas,temperature,pressure);

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

  // World
  //
  G4Box*
  sWorld = new G4Box("World",					//name
                   worldSizeX/2,worldSizeYZ/2,worldSizeYZ/2);	//dimensions

  G4LogicalVolume*		   			                      
  lWorld = new G4LogicalVolume(sWorld,			//shape
                               worldMaterial,		//material
                              "World");			//name

  G4VPhysicalVolume*                                   
  pWorld = new G4PVPlacement(0,				//no rotation
  			   G4ThreeVector(),		//at (0,0,0)
                           lWorld,			//logical volume
			   "World",			//name
                           0,	       		        //mother  volume
                           false,			//no boolean operation
                           0);				//copy number
  //			   
  // Absorber
  //			   
  G4Box*
  sAbsor = new G4Box("Absorber",				//name
                   absorSizeX/2,absorSizeYZ/2,absorSizeYZ/2);	//dimensions
		   			                      
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
  //
  // Tallies (optional)
  //
  if (tallyNumber > 0) {      
    G4Box*
    sTally = new G4Box("Tally",tallySize.x()/2,tallySize.y()/2,tallySize.z()/2);
    lTally = new G4LogicalVolume(sTally,tallyMaterial,"Tally");
    
    for (G4int j=0; j<tallyNumber; j++)
       {
        new G4PVPlacement(0,				//no rotation
  			  tallyPosition[j],		//position
                          lTally,			//logical volume
			  "Tally",			//name
                          lAbsor,	       		//mother  volume
                          false,			//no boolean operation
                          j);				//copy number
       }
       
    tallyMass = tallySize.x()*tallySize.y()*tallySize.z()
               *(tallyMaterial->GetDensity());
  } 

  PrintParameters();
    
  //
  //always return the World volume
  //  
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(absorSizeX,"Length")
         << " of " << absorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
  
  if (tallyNumber > 0) {
    G4cout << "---> There are " << tallyNumber << " tallies : "
           << G4BestUnit(tallySize,"Length")
	   << " of " << tallyMaterial->GetName()
	   << "  (mass : " << G4BestUnit(tallyMass,"Mass") << ")" << G4endl;
	   
    for (G4int j=0; j<tallyNumber; j++)
     G4cout << "tally " << j << ": "
            << "position = " << G4BestUnit(tallyPosition[j],"Length") << G4endl;
    G4cout << "\n---------------------------------------------------------\n";
  }	    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double value)
{
  absorSizeX = value; worldSizeX = 1.2*absorSizeX;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeYZ(G4double value)
{
  absorSizeYZ = value; 
  worldSizeYZ = 1.2*absorSizeYZ;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial) {
    absorMaterial = pttoMaterial;
    if(lAbsor) {
      lAbsor->SetMaterial(absorMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if (magField) delete magField;	//delete the existing magn field
  
  if (fieldValue!=0.)			// create a new one if non nul
    {
      magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));        
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

void DetectorConstruction::SetTallySize(G4ThreeVector value)
{
  tallySize = value;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallyMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (pttoMaterial) {
    tallyMaterial = pttoMaterial;
    if(lTally) {
      lTally->SetMaterial(tallyMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTallyPosition(G4ThreeVector value)
{
  if (tallyNumber < MaxTally) {
    tallyPosition[tallyNumber] = value; 
    tallyNumber++;
  }
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
