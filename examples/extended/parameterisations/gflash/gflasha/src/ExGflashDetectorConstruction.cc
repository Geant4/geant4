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
//
/// \file ExGflashDetectorConstruction.cc
/// \brief Implementation of the ExGflashDetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// User Classes
#include "ExGflashDetectorConstruction.hh"
#include "ExGflashSensitiveDetector.hh"
#include "ExGflashMessenger.hh"
#include "ExGflashHomoShowerTuning.hh"

// G4 Classes
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"
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

G4ThreadLocal GFlashShowerModel* ExGflashDetectorConstruction::fFastShowerModel = 0; 
G4ThreadLocal 
GFlashHomoShowerParameterisation* ExGflashDetectorConstruction::fParameterisation = 0;
G4ThreadLocal GFlashParticleBounds* ExGflashDetectorConstruction::fParticleBounds = 0;
G4ThreadLocal GFlashHitMaker* ExGflashDetectorConstruction::fHitMaker = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int kMaxBin = 500;
//const G4int kMaxRowCol = 50;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashDetectorConstruction::ExGflashDetectorConstruction()
  :G4VUserDetectorConstruction(),fNbOfCrystals(10),
   fCrystal_log(0),fDetMat(0),fRegion(0),
   fVerbose(0),fNLtot(40),fNRtot(50),fDLradl(0.5),fDRradl(0.1)
{
  G4cout<<"ExGflashDetectorConstruction::Detector constructor"<<G4endl;
  fGflashMessenger = new ExGflashMessenger(this);

  // Crystall
  fCrystalWidth = 3*cm;
  fCrystalLength = 24*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashDetectorConstruction::~ExGflashDetectorConstruction()
{
  delete fGflashMessenger;
  delete fFastShowerModel; 
  delete fParameterisation;
  delete fParticleBounds;
  delete fHitMaker;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExGflashDetectorConstruction::Construct()
{
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  G4cout << "Defining the materials" << G4endl;
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  // Build materials
  G4Material* air   = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* pbWO4 = nistManager->FindOrBuildMaterial("G4_PbWO4");

  fDetMat = pbWO4;
  
  //------------------------------ 
  // Calorimeter segments
  //------------------------------
  // Simplified `CMS-like` PbWO4 crystal calorimeter  
  
                       
  // Simplified `CMS-like` PbWO4 crystal calorimeter  

  // Calorimeter 
  G4double calo_xside = (fCrystalWidth*fNbOfCrystals);
  G4double calo_yside = (fCrystalWidth*fNbOfCrystals);
  G4double calo_zside = fCrystalLength;

  //The Experimental Hall
  G4double experimentalHall_x=calo_xside * 4;
  G4double experimentalHall_y=calo_yside * 4;
  G4double experimentalHall_z=calo_zside * 4;
  
  G4VSolid* experimentalHall_box 
    = new G4Box("expHall_box",             // World Volume
                experimentalHall_x,        // x size
                experimentalHall_y,        // y size
                experimentalHall_z);       // z size
  
  G4LogicalVolume* experimentalHall_log 
    = new G4LogicalVolume(experimentalHall_box,
                          air,
                          "expHall_log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  G4VPhysicalVolume* experimentalHall_phys 
    = new G4PVPlacement(0,
                        G4ThreeVector(),   //at (0,0,0)
                        "expHall",
                        experimentalHall_log,
                        0,
                        false,
                        0);
  
  G4Box* calo_box= new G4Box("Calorimeter",  // its name
                             calo_xside/2.,      // size
                             calo_yside/2.,
                             calo_zside/2.);
  G4LogicalVolume* calo_log 
    = new G4LogicalVolume(calo_box,      // its solid
                          air,           // its material
                          "calo_log",    // its name
                          0,             // opt: fieldManager
                          0,             // opt: SensitiveDetector 
                          0);            // opt: UserLimit
  
  G4double xpos = 0.0;
  G4double ypos = 0.0;
  G4double zpos = calo_zside/2.;// face @ z= 0.0

  new G4PVPlacement(0,
                    G4ThreeVector(xpos, ypos, zpos),
                    calo_log,
                    "calorimeter",
                    experimentalHall_log,
                    false,
                    1);
          
  // Crystals
  G4VSolid* crystal_box 
    = new G4Box("Crystal",                // its name
                 fCrystalWidth/2,
                 fCrystalWidth/2,
                 fCrystalLength/2);   
                           // size
  fCrystal_log 
    = new G4LogicalVolume(crystal_box,    // its solid
                          fDetMat,          // its material
                          "Crystal_log"); // its name
          
  for (G4int i=0; i<fNbOfCrystals; i++)
    {
      
      for (G4int j=0; j<fNbOfCrystals; j++)
        {  
          G4int n =  i*10+j;
          G4ThreeVector crystalPos((i*fCrystalWidth)-(calo_xside-fCrystalWidth)/2.,
                                   (j*fCrystalWidth)-(calo_yside-fCrystalWidth)/2.,0 );
            new G4PVPlacement(0,               // no rotation
                                crystalPos,      // translation
                                fCrystal_log,
                                "crystal",       // its name
                                calo_log,
                                false,
                                n);
        }
    }  
  G4cout << "There are " << fNbOfCrystals <<
    " crystals per row in the calorimeter, so in total "<<
    fNbOfCrystals*fNbOfCrystals << " crystals" << G4endl;
  G4cout << "Total Calorimeter size" << calo_xside / cm << " cm x "
         << calo_yside / cm << " cm x "
         << calo_zside / cm << " cm"  << G4endl;
  G4cout << "They have width of  " << fCrystalWidth /cm <<
    "  cm and a length of  " <<  fCrystalLength /cm
         <<" cm. The Material is "<< fDetMat 
         <<" Total: "<< fCrystalLength / fDetMat->GetRadlen()
         <<" X0" << G4endl;
  
  
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* caloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  calo_log->SetVisAttributes(caloVisAtt);
  fCrystal_log->SetVisAttributes(crystalVisAtt);

  // define the fParameterisation region
  G4cout << "\n ---> DetectorConstruction Region Definition" << G4endl;
  fRegion = new G4Region("crystals");
  calo_log->SetRegion(fRegion);
  fRegion->AddRootLogicalVolume(calo_log);
  
    return experimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashDetectorConstruction::ConstructSDandField()
{
  // -- sensitive detectors:

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  ExGflashSensitiveDetector* SD =
    new ExGflashSensitiveDetector("Calorimeter",this);

  SDman->AddNewDetector(SD);
  if (fCrystal_log)
    fCrystal_log->SetSensitiveDetector(SD);
  
  // -- fast simulation models:
  // **********************************************
  // * Initializing shower modell
  // ***********************************************
  G4cout << "\n--> Creating shower parameterization models" << G4endl;
  fFastShowerModel = new GFlashShowerModel("fFastShowerModel", fRegion);
  fParameterisation =
    new GFlashHomoShowerParameterisation(fDetMat,
                                         new ExGflashHomoShowerTuning());
  fFastShowerModel->SetParameterisation(*fParameterisation);
  // Energy Cuts to kill particles:
  fParticleBounds = new GFlashParticleBounds();
  fFastShowerModel->SetParticleBounds(*fParticleBounds);
  // Makes the EnergieSpots
  fHitMaker = new GFlashHitMaker();
  fFastShowerModel->SetHitMaker(*fHitMaker);
  G4cout<<"end shower parameterization."<<G4endl;
  // **********************************************
  // Get Rad Len and R molere ?
  fSDRadLen = fDetMat->GetRadlen();
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashDetectorConstruction::SetLBining(G4ThreeVector Value)
{
  fNLtot = (G4int)Value(0);
  if (fNLtot > kMaxBin) {
    G4cout << "\n ---> warning from SetLBining: "
           << fNLtot << " truncated to " << kMaxBin << G4endl;
    fNLtot = kMaxBin;
  }  
  fDLradl = Value(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashDetectorConstruction::SetRBining(G4ThreeVector Value)
{
  fNRtot = (G4int)Value(0);
  if (fNRtot > kMaxBin) {
    G4cout << "\n ---> warning from SetRBining: "
           << fNRtot << " truncated to " << kMaxBin << G4endl;
    fNRtot = kMaxBin;
  }    
  fDRradl = Value(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashDetectorConstruction::SetMaterial(G4String mat)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if(pttoMaterial &&  fDetMat != pttoMaterial) {
    fDetMat = pttoMaterial;
    if(fCrystal_log) { fCrystal_log->SetMaterial(fDetMat); }
    if(fParameterisation) {
      fParameterisation->SetMaterial(fDetMat);
      fParameterisation->PrintMaterial(fDetMat);
      // Get Rad Len and R molere ?
      fSDRadLen = fDetMat->GetRadlen();
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
