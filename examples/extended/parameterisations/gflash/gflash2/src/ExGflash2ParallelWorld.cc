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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// User Classes
#include "ExGflash2ParallelWorld.hh"

// G4 Classes
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoDelete.hh"
#include "globals.hh"

//fast simulation
#include "GFlashHomoShowerParameterisation.hh"
#include "G4FastSimulationManager.hh"
#include "GFlashShowerModel.hh"
#include "GFlashHitMaker.hh"
#include "GFlashParticleBounds.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal GFlashShowerModel* ExGflash2ParallelWorld::fFastShowerModel = 0; 
G4ThreadLocal 
GFlashHomoShowerParameterisation* ExGflash2ParallelWorld::fParameterisation = 0;
G4ThreadLocal GFlashParticleBounds* ExGflash2ParallelWorld::fParticleBounds = 0;
G4ThreadLocal GFlashHitMaker* ExGflash2ParallelWorld::fHitMaker = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflash2ParallelWorld::ExGflash2ParallelWorld(G4String aWorldName)
  :G4VUserParallelWorld(aWorldName), fRegion(nullptr)
{
  G4cout<<"ExGflash2ParallelWorld::Parralel world constructor"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflash2ParallelWorld::~ExGflash2ParallelWorld()
{
  delete fFastShowerModel; 
  delete fParameterisation;
  delete fParticleBounds;
  delete fHitMaker;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflash2ParallelWorld::Construct()
{
  // In parallel world material does not matter
  G4Material* dummy   = nullptr;

  // Build parallel/ghost geometry:
  auto ghostLogicalVolume = GetWorld()->GetLogicalVolume();
  
  // Use part of the Ex1GflashDetectorConstruction (without individual crystals)
  
  G4int nbOfCrystals = 10;  // this are the crystals PER ROW in this example 
                       // cube of 10 x 10 crystals 
                       // don't change it @the moment, since 
                       // the readout in event action assumes this 
                       // dimensions and is not automatically adapted
                       // in this version of the example :-( 
  // Simplified `CMS-like` PbWO4 crystal calorimeter  
  G4double calo_xside = 31*cm;
  G4double calo_yside = 31*cm;
  G4double calo_zside = 24*cm; 

  G4double crystalWidth = 3*cm;
  G4double crystalLength = 24*cm;

  calo_xside = (crystalWidth*nbOfCrystals)+1*cm;
  calo_yside = (crystalWidth*nbOfCrystals)+1*cm;
  calo_zside = crystalLength;
  
  G4Box* calo_box= new G4Box("CMS calorimeter",  // its name
                             calo_xside/2.,      // size
                             calo_yside/2.,
                             calo_zside/2.);
  G4LogicalVolume* calo_log 
    = new G4LogicalVolume(calo_box,      // its solid
                          dummy,           // its material
                          "calo log",    // its name
                          0,             // opt: fieldManager
                          0,             // opt: SensitiveDetector 
                          0);            // opt: UserLimit

  G4double xpos = 0.0;
  G4double ypos = 0.0;
  G4double zpos = 100.0*cm;
  new G4PVPlacement(0,
                    G4ThreeVector(xpos, ypos, zpos),
                    calo_log,
                    "calorimeter",
                    ghostLogicalVolume,
                    false,
                    1);

  G4VisAttributes* caloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  calo_log->SetVisAttributes(caloVisAtt);

  // define the fParameterisation region
  fRegion = new G4Region("crystals");
  calo_log->SetRegion(fRegion);
  fRegion->AddRootLogicalVolume(calo_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflash2ParallelWorld::ConstructSD()
{
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material*          pbWO4 = nistManager->FindOrBuildMaterial("G4_PbWO4");
  // -- fast simulation models:
  // **********************************************
  // * Initializing shower modell
  // ***********************************************
  G4cout << "Creating shower parameterization models" << G4endl;
  fFastShowerModel = new GFlashShowerModel("fFastShowerModel", fRegion);
  fParameterisation = new GFlashHomoShowerParameterisation(pbWO4);
  fFastShowerModel->SetParameterisation(*fParameterisation);
  // Energy Cuts to kill particles:
  fParticleBounds = new GFlashParticleBounds();
  fFastShowerModel->SetParticleBounds(*fParticleBounds);
  // Makes the EnergieSpots
  fHitMaker = new GFlashHitMaker();
  fFastShowerModel->SetHitMaker(*fHitMaker);
  G4cout<<"end shower parameterization."<<G4endl;
  // **********************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
