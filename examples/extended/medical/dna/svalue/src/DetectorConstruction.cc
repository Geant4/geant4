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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 78723 2014-01-20 10:32:17Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  //default tracking cut  
  fTrackingCut = 7.4*CLHEP::eV;
  
  // default parameter values
  fWorldRadius = 10*CLHEP::m;
  fCytoThickness = 20*CLHEP::nm;
  fNuclRadius = 10*CLHEP::nm;
  
  DefineMaterials();

  // create commands for interactive definition of the detector  
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
  G4NistManager* man = G4NistManager::Instance();
  
  fWorldMaterial = man->FindOrBuildMaterial("G4_WATER");
  fCytoMaterial = fNuclMaterial = man->FindOrBuildMaterial("G4_WATER");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if(fWorld) return fWorld;
                   
  // Spherical world
  //
  G4Sphere* 
  sWorld = new G4Sphere("World",                           
                        0., 
                        1000*fNuclRadius, 
                        0., 
                        twopi, 
                        0., 
                        pi);    

  fLogicalWorld = new G4LogicalVolume(sWorld,                        
                                      fWorldMaterial,          
                                      "World");              
                                   
  fWorld = new G4PVPlacement(0,                         
                             G4ThreeVector(),           
                             fLogicalWorld,                    
                             "World",                 
                             0,                         
                             false,                     
                             0);                        

  // Spherical nucleus
  //
  G4Sphere* 
  sNucl = new G4Sphere("Nucl",                           
                       0., 
                       fNuclRadius, 
                       0., 
                       twopi, 
                       0., 
                       pi);    

  fLogicalNucl = new G4LogicalVolume(sNucl,                        
                                     fNuclMaterial,          
                                    "Nucl");              
                                   
  fNucl = new G4PVPlacement(0,                         
                            G4ThreeVector(),          
                            "Nucl",                 
                            fLogicalNucl,                         
                            fWorld,
                            false,                     
                            0);                       

  // Spherical shell for cytoplasm
  //
  G4Sphere* 
  sCyto = new G4Sphere("Cyto",                           
                        fNuclRadius, 
                        fNuclRadius+fCytoThickness, 
                        0., 
                        twopi, 
                        0., 
                        pi);

  fLogicalCyto = new G4LogicalVolume(sCyto,                    
                                     fCytoMaterial,       
                                    "Cyto");       
                                   
  fCyto = new G4PVPlacement(0,                         
                            G4ThreeVector(),           
                            "Cyto",                
                            fLogicalCyto,
                            fWorld,                       
                            false,                   
                            0);               
  //
  
  PrintParameters();
    
  // Tracking cut
  //
  fLogicalNucl->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));
  fLogicalCyto->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));
  fLogicalWorld->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));
    
  //
  //always return the root volume
  //  
  return fWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() const
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The tracking cut is set to " 
         << G4BestUnit(fTrackingCut,"Energy") << G4endl;
  G4cout << "---> The World is a sphere of " 
         << G4BestUnit(1000*fNuclRadius,"Length") << "radius of "
         << fWorldMaterial->GetName() << G4endl;
  G4cout << "---> The Nucleus is a sphere of " 
         << G4BestUnit(fNuclRadius,"Length") << "radius of "
         << fWorldMaterial->GetName() << " of mass "
  << G4BestUnit(GetNuclMass(),"Mass") << G4endl;
  G4cout << "---> The Cytoplasm is a spherical shell of thickness " 
         << G4BestUnit(fCytoThickness,"Length") << "of "
         << fWorldMaterial->GetName() << " of mass "
  << G4BestUnit(GetCytoMass(),"Mass") << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTrackingCut(G4double value)
{
  fTrackingCut = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNuclRadius(G4double value)
{
  fNuclRadius = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCytoThickness(G4double value)
{
  fCytoThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && pttoMaterial != fWorldMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLogicalWorld) fLogicalWorld->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCytoMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && pttoMaterial != fCytoMaterial) {
    fCytoMaterial = pttoMaterial;
    if(fLogicalCyto) fLogicalCyto->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNuclMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && pttoMaterial != fNuclMaterial) {
    fNuclMaterial = pttoMaterial;
    if(fLogicalNucl) fLogicalNucl->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
