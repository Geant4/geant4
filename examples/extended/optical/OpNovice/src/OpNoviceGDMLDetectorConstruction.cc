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

#include "OpNoviceGDMLDetectorConstruction.hh"
#include "OpNoviceDetectorMessenger.hh"

#include "globals.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OpNoviceGDMLDetectorConstruction::OpNoviceGDMLDetectorConstruction(
  G4String fname)
  : G4VUserDetectorConstruction()
{
  fDumpGdmlFileName = "OpNovice_dump.gdml";
  fVerbose          = false;
  fDumpGdml         = false;
  fGdmlFile         = fname;
  // create a messenger for this class
  fDetectorMessenger = new OpNoviceDetectorMessenger(this);

  G4cout << "Building detector from GDML file: " << fname << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OpNoviceGDMLDetectorConstruction::~OpNoviceGDMLDetectorConstruction()
{
  delete fDetectorMessenger;
  delete fParser;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* OpNoviceGDMLDetectorConstruction::Construct()
{
  ReadGDML();
  G4VPhysicalVolume* worldPhysVol = fParser->GetWorldVolume();
  if(fDumpGdml)
    fParser->Write(fDumpGdmlFileName, worldPhysVol);
  return worldPhysVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::ConstructSDandField() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::ReadGDML()
{
  fParser = new G4GDMLParser();
  fParser->Read(fGdmlFile, false);
  G4VPhysicalVolume* world = fParser->GetWorldVolume();
  // GDML parser makes world invisible. make it visible again.
  G4LogicalVolume* pworldLogical = world->GetLogicalVolume();
  pworldLogical->SetVisAttributes(0);
  G4cout << world->GetTranslation() << G4endl << G4endl;
  if(fVerbose)
  {
    G4cout << "Found world:  " << world->GetName() << G4endl;
    G4cout << "world LV:  " << world->GetLogicalVolume()->GetName() << G4endl;
  }
  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  if(fVerbose)
  {
    G4cout << "Found " << pLVStore->size() << " logical volumes." << G4endl
           << G4endl;
  }
  G4PhysicalVolumeStore* pPVStore = G4PhysicalVolumeStore::GetInstance();
  if(fVerbose)
  {
    G4cout << "Found " << pPVStore->size() << " physical volumes." << G4endl
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::SetDumpGdml(G4bool val)
{
  fDumpGdml = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool OpNoviceGDMLDetectorConstruction::IsDumpGdml() const
{
  return fDumpGdml;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::SetVerbose(G4bool val)
{
  fVerbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool OpNoviceGDMLDetectorConstruction::IsVerbose() const { return fVerbose; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::SetDumpGdmlFile(G4String val)
{
  fDumpGdmlFileName = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String OpNoviceGDMLDetectorConstruction::GetDumpGdmlFileName() const
{
  return fDumpGdmlFileName;
}
