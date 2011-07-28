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

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4RunManager.hh" 
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(RunAction* r, PrimaryGeneratorAction* p)
  : run(r), gen(p)
{
  physWorld = 0;
  physAbs   = 0;
  SetWorldMaterial   ("G4_Galactic");
  SetAbsorberMaterial("G4_Fe");
  SetAbsorberWidth(0.24*mm);
  SetBeamMomentum(172.*MeV);
  SetBeamMomentumSpread(1.*MeV);
  SetFileName("test41.log");
  
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

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  // Cleanup old geometry
  if(physWorld) {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
  }
  // World
  G4double y = 10*cm;
  G4Box* solidWorld = new G4Box("World",y + 1.,y + 1.,width + 1);	
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,matWorld,"World");	
  physWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0); 

  // Absorber
  G4Box* solidAbs = new G4Box("Abs",y,y,0.5*width);	
  G4LogicalVolume* logicAbs = new G4LogicalVolume(solidAbs,matAbsorber,"Abs");	
  physAbs = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.5*width),
			      logicAbs,"Abs",logicWorld,false,0); 

  G4cout << "### Absorber width(mm)= " << width/mm 
	 << " made of " << matAbsorber->GetName()
	 << G4endl;

  return physWorld;
}

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* m = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(m) {
    matWorld = m;
    if(physWorld) G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

void DetectorConstruction::SetAbsorberMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* m = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(m) {
    matAbsorber = m;
    run->SetMaterialName(mat);
    if(physWorld) G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

void DetectorConstruction::SetAbsorberWidth(G4double d)
{
  width = d;
  run->SetAbsorberWidth(width);
  if(physWorld) G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
 
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
