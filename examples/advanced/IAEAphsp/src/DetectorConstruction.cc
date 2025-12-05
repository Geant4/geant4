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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  G4NistManager* man = G4NistManager::Instance();
  //man->SetVerbose(1);

  fMessenger = new DetectorMessenger(this);

  fWorldMat   = man->FindOrBuildMaterial("G4_Galactic");
  fWorldXY = 50.*cm;
  fWorldZ = 100.*cm;

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DumpGeometryParameters();

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  //
  // World
  //

  G4Box* sWorld = new G4Box("World",
			    fWorldXY, fWorldXY, fWorldZ);
  G4LogicalVolume* lWorld = new G4LogicalVolume(sWorld,
						fWorldMat, "World");
  G4VPhysicalVolume* phWorld = new G4PVPlacement(0, G4ThreeVector(),
						 "World", lWorld,
						 0, false, 0);

  //
  // Visualization attributes
  //
  auto visAtt = new G4VisAttributes(G4Colour(0.1,0.5,1.0));
  visAtt->SetVisibility(true);
  lWorld->SetVisAttributes(visAtt);

  return phWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DumpGeometryParameters()
{

  G4cout << "\n===================================================" << G4endl;
  G4cout << "#               IAEAphsp Geometry                 #" << G4endl;
  G4cout << "===================================================" << G4endl;
  G4cout << "  WorldXY = " << fWorldXY/cm << " cm " << G4endl;
  G4cout << "  WorldZ = " << fWorldZ/cm << " cm " << G4endl;
  G4cout << "  WorldMat: " << fWorldMat->GetName() << G4endl;
  G4cout << "===================================================\n" << G4endl;
}
