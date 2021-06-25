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
// Geant4 headers
#include "G4RunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4GDMLParser.hh"
// project headers
#include "OpNoviceGDMLDetectorConstruction.hh"
#include "OpNoviceGDMLDetectorConstructionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OpNoviceGDMLDetectorConstruction::OpNoviceGDMLDetectorConstruction(G4String fname)
: G4VUserDetectorConstruction() {
    fDumpgdmlFile = "OpNovice_dump.gdml";
    fverbose = false;
    fdumpgdml = false;
    gdmlFile = fname;
    // create a messenger for this class
    fDetectorMessenger = new OpNoviceGDMLDetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OpNoviceGDMLDetectorConstruction::~OpNoviceGDMLDetectorConstruction() {
    delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume *OpNoviceGDMLDetectorConstruction::Construct() {
    ReadGDML();
    G4VPhysicalVolume *worldPhysVol = parser->GetWorldVolume();
    if (fdumpgdml) parser->Write(fDumpgdmlFile, worldPhysVol);
    return worldPhysVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::ConstructSDandField() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::ReadGDML() {
    parser = new G4GDMLParser();
    parser->Read(gdmlFile, false);
    G4VPhysicalVolume *world = parser->GetWorldVolume();
    //----- GDML parser makes world invisible, this is a hack to make it
    //visible again...
    G4LogicalVolume *pworldLogical = world->GetLogicalVolume();
    pworldLogical->SetVisAttributes(0);
    G4cout << world->GetTranslation() << G4endl << G4endl;
    if (fverbose) {
        G4cout << "Found world:  " << world-> GetName() << G4endl;
        G4cout << "world LV:  " << world->GetLogicalVolume()->GetName() << G4endl;
    }
    G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
    if (fverbose) {
        G4cout << "Found " << pLVStore->size()
                << " logical volumes."
                << G4endl << G4endl;
    }
    G4PhysicalVolumeStore *pPVStore = G4PhysicalVolumeStore::GetInstance();
    if (fverbose) {
        G4cout << "Found " << pPVStore->size()
                << " physical volumes."
                << G4endl << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::UpdateGeometry() {
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstruction::SetDumpgdml(G4bool fdumpgdml1) {
    this->fdumpgdml = fdumpgdml1;
}
G4bool OpNoviceGDMLDetectorConstruction::IsDumpgdml() const {
    return fdumpgdml;
}
void OpNoviceGDMLDetectorConstruction::SetVerbose(G4bool fverbose1) {
    this->fverbose = fverbose1;
}
G4bool OpNoviceGDMLDetectorConstruction::IsVerbose() const {
    return fverbose;
}
void OpNoviceGDMLDetectorConstruction::SetDumpgdmlFile(G4String fDumpgdmlFile1) {
    this->fDumpgdmlFile = fDumpgdmlFile1;
}
G4String OpNoviceGDMLDetectorConstruction::GetDumpgdmlFile() const {
    return fDumpgdmlFile;
}