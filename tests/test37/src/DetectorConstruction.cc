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
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4RunManager.hh" 

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:solidWorld(0),logicWorld(0),physiWorld(0),
 solidMedium1(0),logicMedium1(0),physiMedium1(0),
 solidLayerMedium1(0),logicLayerMedium1(0),physiLayerMedium1(0),
 solidMedium2(0),logicMedium2(0),physiMedium2(0),
 solidLayerMedium2(0),logicLayerMedium2(0),physiLayerMedium2(0),
 solidLayerMedium3(0),logicLayerMedium3(0),physiLayerMedium3(0)
{
  // default parameter values of the calorimeter
  Absorber1Thickness = 0.500*(0.325/16.6)*cm;// 0.500 FMR
  Absorber2Thickness = 0.747*(0.569/2.7)*cm; // 0.747 FMR
  Absorber3Thickness = 0.747*(0.569/2.7)*cm; // 0.747 FMR
  NbOfLayersOfMedium1        = 1;
  NbOfLayersOfMedium2        = 1;
  NbOfLayersOfMedium3        = 1;
  LayerThichness1 = Absorber1Thickness/NbOfLayersOfMedium1; 
  LayerThichness2 = Absorber2Thickness/NbOfLayersOfMedium2;
  LayerThichness3 = Absorber3Thickness/NbOfLayersOfMedium3;

  // materials  
  DefineMaterials();
  SetWorldMaterial   ("G4_Galactic");
  SetAbsorber1Material("G4_Galactic");
  SetAbsorber2Material("G4_Galactic");
  SetAbsorber3Material("G4_Galactic");
  
  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new DetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterials()
{ 
  // This function illustrates the possible ways to define materials using 
  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(1);

}
G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  DefineMaterials();
  // World
  solidWorld = new G4Box("World",(4./2)*cm,(15./2)*cm,(15./2)*cm);	
  logicWorld = new G4LogicalVolume(solidWorld,WorldMaterial,"World");	
  physiWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0);			

//------------------ Definition of Mediums 1 & 2 & 3 ------------------------------
//--------------------------------------------------------------------------------
//                                 Medium1

  solidMedium1 = new G4Box("Medium 1",Absorber1Thickness/2,(10.0/2)*cm,(10./2)*cm);
  logicMedium1 = new G4LogicalVolume(solidMedium1,Absorber1Material,"Medium 1");	
  physiMedium1 = new G4PVPlacement(0,G4ThreeVector(Absorber1Thickness/2,0.,0.),logicMedium1,
                                   "Medium 1",logicWorld,false,0);			

                                // Layers of Medium1
  LayerThichness1 = Absorber1Thickness/NbOfLayersOfMedium1;
  solidLayerMedium1 = new G4Box("LayerMedium 1",LayerThichness1/2,(10.0/2)*cm,(10./2)*cm);
  logicLayerMedium1 = new G4LogicalVolume(solidLayerMedium1,Absorber1Material,"LayerMedium 1");	
  
  for (G4int iAbs=0; iAbs<NbOfLayersOfMedium1; iAbs++) 
  {
    physiLayerMedium1 = new G4PVPlacement(0,
	G4ThreeVector(((LayerThichness1-Absorber1Thickness)/2 + iAbs * LayerThichness1)  
		      ,0.,0.),logicLayerMedium1,"LayerMedium 1",logicMedium1,false,iAbs+1);			
  }
//--------------------------------------------------------------------------------
//                                 Medium2
  solidMedium2 = new G4Box("Medium 2",Absorber2Thickness/2,(10.0/2)*cm,(10./2)*cm);
  logicMedium2 = new G4LogicalVolume(solidMedium2,Absorber2Material,"Medium 2");	
  physiMedium2 = new G4PVPlacement(0,
               G4ThreeVector(Absorber1Thickness+0.5*Absorber2Thickness,0.,0.),
	       logicMedium2,"Medium 2",logicWorld,false,0);			

                                // Layers of Medium2
  LayerThichness2 = Absorber2Thickness/NbOfLayersOfMedium2;
  solidLayerMedium2 = new G4Box("LayerMedium 2",LayerThichness2/2,(10.0/2)*cm,(10./2)*cm);
  logicLayerMedium2 = new G4LogicalVolume(solidLayerMedium2,Absorber2Material,"LayerMedium 2");	
  
  for (G4int iAb=0; iAb<NbOfLayersOfMedium2; iAb++) 
  {
    physiLayerMedium2 = new G4PVPlacement(0,
	      G4ThreeVector(((LayerThichness2-Absorber2Thickness)/2 + iAb * LayerThichness2)  
			    ,0.,0.),logicLayerMedium2,"LayerMedium 2",logicMedium2,false,iAb+1);			
  }
//--------------------------------------------------------------------------------
//                                 Medium3
  solidMedium3 = new G4Box("Medium 3",Absorber3Thickness/2,(10.0/2)*cm,(10./2)*cm);
  logicMedium3 = new G4LogicalVolume(solidMedium3,Absorber3Material,"Medium 3");	
  physiMedium3 = new G4PVPlacement(0,
               G4ThreeVector(Absorber1Thickness+Absorber2Thickness+0.5*Absorber3Thickness,0.,0.),
	       logicMedium3,"Medium 3",logicWorld,false,0);			

                                // Layers of Medium2
  LayerThichness3 = Absorber3Thickness/NbOfLayersOfMedium3;
  solidLayerMedium3 = new G4Box("LayerMedium 3",LayerThichness3/2,(1.0/2)*cm,(10./2)*cm);
  logicLayerMedium3 = new G4LogicalVolume(solidLayerMedium3,Absorber3Material,"LayerMedium 3");	
  
  for (G4int iAb=0; iAb<NbOfLayersOfMedium3; iAb++) 
  {
    physiLayerMedium3 = new G4PVPlacement(0,
               G4ThreeVector(((LayerThichness3-Absorber3Thickness)/2 + iAb * LayerThichness3)  
			     ,0.,0.),logicLayerMedium3,"LayerMedium 3",logicMedium3,false,iAb+1);			
  }

  return physiWorld;
}

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);     
  if (pttoMaterial) WorldMaterial = pttoMaterial;
}


void DetectorConstruction::SetAbsorber1Material(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);     
  if (pttoMaterial) Absorber1Material = pttoMaterial;                  
}


void DetectorConstruction::SetAbsorber2Material(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);     
  if (pttoMaterial) Absorber2Material = pttoMaterial;                  
}

void DetectorConstruction::SetAbsorber3Material(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);     
  if (pttoMaterial) Absorber3Material = pttoMaterial;                  
}

void DetectorConstruction::SetAbsorber1Thickness(G4double val)
{
  Absorber1Thickness = val;
}  

void DetectorConstruction::SetAbsorber2Thickness(G4double val)
{
  Absorber2Thickness = val;
}  

void DetectorConstruction::SetAbsorber3Thickness(G4double val)
{
  Absorber3Thickness = val;
}  

void DetectorConstruction::SetNbOfLayersOfMedium1(G4int val)
{
  NbOfLayersOfMedium1 = val;
}

void DetectorConstruction::SetNbOfLayersOfMedium2(G4int val)
{
  NbOfLayersOfMedium2 = val;
}

void DetectorConstruction::SetNbOfLayersOfMedium3(G4int val)
{
  NbOfLayersOfMedium3 = val;
}
 
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}


