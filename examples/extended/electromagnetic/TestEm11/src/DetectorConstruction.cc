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
// $Id: DetectorConstruction.cc,v 1.9 2009-11-16 15:09:10 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  // default parameter values
  absorSizeX = absorSizeYZ = 1*mm;
  nbOfLayers = 1;
  
  absorMaterial = 0;
  magField = 0;
  pAbsor   = 0;
  
  DefineMaterials();
  SetMaterial("G4_Si");

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
  G4NistManager* man = G4NistManager::Instance();
  
  G4bool isotopes = false;
  
  man->FindOrBuildMaterial("G4_Al", isotopes);
  man->FindOrBuildMaterial("G4_Si", isotopes);
  man->FindOrBuildMaterial("G4_Fe", isotopes);
  man->FindOrBuildMaterial("G4_Cu", isotopes);  
  man->FindOrBuildMaterial("G4_Ge", isotopes);
  man->FindOrBuildMaterial("G4_Mo", isotopes);
  man->FindOrBuildMaterial("G4_Ta", isotopes);
  man->FindOrBuildMaterial("G4_W" , isotopes);
  man->FindOrBuildMaterial("G4_Au", isotopes);
  man->FindOrBuildMaterial("G4_Pb", isotopes);  
  man->FindOrBuildMaterial("G4_PbWO4", isotopes);
  man->FindOrBuildMaterial("G4_SODIUM_IODIDE", isotopes);
  
  man->FindOrBuildMaterial("G4_AIR"  , isotopes);
  ///man->FindOrBuildMaterial("G4_WATER", isotopes);
  
  G4Element* H = man->FindOrBuildElement("H", isotopes); 
  G4Element* O = man->FindOrBuildElement("O", isotopes);
  
  G4Material* H2O = 
  new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
        
 G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
		   
  // Absorber
  //
  G4Box* 
  sAbsor = new G4Box("Absorber",				//name
    		 absorSizeX/2,absorSizeYZ/2,absorSizeYZ/2);	//size

  G4LogicalVolume*
  lAbsor = new G4LogicalVolume(sAbsor,			//solid
      			       absorMaterial,		//material
      			      "Absorber");		//name
				   
  pAbsor = new G4PVPlacement(0,				//no rotation
                             G4ThreeVector(),		//at (0,0,0)
                             lAbsor,			//logical volume
                            "Absorber",			//name
                             0,				//mother  volume
                             false,			//no boolean operation
                             0);			//copy number

  // Layers
  //
  layerThickness = absorSizeX/nbOfLayers;
  
  if (nbOfLayers > 1)  {                           
  			   
    G4Box*
    sLayer = new G4Box("Layer",					//name
                 layerThickness/2,absorSizeYZ/2,absorSizeYZ/2);	//dimensions
		 
    G4LogicalVolume*		   			                      
    lLayer = new G4LogicalVolume(sLayer,		//shape
                                 absorMaterial,		//material
                                 "Layer");		//name
			                           
             new G4PVReplica("Layer",			//name
                           lLayer,			//logical volume
                           lAbsor,	       		//mother  volume
                           kXAxis,			//axis   of replica
                           nbOfLayers,			//number of replica
                           layerThickness);		//width  of replica
  }			   

  PrintParameters();
    
  //
  //always return the root volume
  //  
  return pAbsor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(absorSizeX,"Length")
         << " of " << absorMaterial->GetName() << " divided in " << nbOfLayers 
	 << " slices of " << G4BestUnit(layerThickness,"Length") << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double value)
{
  absorSizeX = value;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeYZ(G4double value)
{
  absorSizeYZ = value; 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) absorMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int value)
{
  nbOfLayers = value; 
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

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

#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()
{
G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
