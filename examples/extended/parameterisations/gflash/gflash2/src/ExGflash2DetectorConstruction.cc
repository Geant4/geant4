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
/// \file ExGflash2DetectorConstruction.cc
/// \brief Implementation of the ExGflash2DetectorConstruction class
//
// Created by Joanna Weng 26.11.2004

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// User Classes
#include "ExGflash2DetectorConstruction.hh"
#include "ExGflash2SensitiveDetector.hh"

// G4 Classes
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



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflash2DetectorConstruction::ExGflash2DetectorConstruction()
  :G4VUserDetectorConstruction(), fCrystalLog(nullptr), fCrystalPhys{}
{
  G4cout<<"ExGflash2DetectorConstruction::Detector constructor"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflash2DetectorConstruction::~ExGflash2DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExGflash2DetectorConstruction::Construct()
{
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  G4cout << "Defining the materials" << G4endl;
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  // Build materials
  G4Material* air   = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* pbWO4 = nistManager->FindOrBuildMaterial("G4_PbWO4");
  
  /*******************************
   * The Experimental Hall       *
   *******************************/
  G4double experimentalHall_x=1000.*cm;
  G4double experimentalHall_y=1000.*cm;
  G4double experimentalHall_z=1000.*cm;
  
  G4VSolid* experimentalHall_box 
    = new G4Box("expHall_box",             // World Volume
                experimentalHall_x,        // x size
                experimentalHall_y,        // y size
                experimentalHall_z);       // z size
  
  G4LogicalVolume* experimentalHallLog 
    = new G4LogicalVolume(experimentalHall_box,
                          air,
                          "expHallLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  G4VPhysicalVolume* experimentalHallPhys 
    = new G4PVPlacement(0,
                        G4ThreeVector(),   //at (0,0,0)
                        "expHall",
                        experimentalHallLog,
                        0,
                        false,
                        0);
  
  
  //------------------------------ 
  // Calorimeter segments
  //------------------------------
  // Simplified `CMS-like` PbWO4 crystal calorimeter  
  
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
  G4LogicalVolume* caloLog 
    = new G4LogicalVolume(calo_box,      // its solid
                          air,           // its material
                          "calo log",    // its name
                          0,             // opt: fieldManager
                          0,             // opt: SensitiveDetector 
                          0);            // opt: UserLimit
  
  G4double xpos = 0.0;
  G4double ypos = 0.0;
  G4double zpos = 100.0*cm;
  new G4PVPlacement(0,
                    G4ThreeVector(xpos, ypos, zpos),
                    caloLog,
                    "calorimeter",
                    experimentalHallLog,
                    false,
                    1);
          
  // Crystals
  G4VSolid* crystal_box 
    = new G4Box("Crystal",                // its name
                 crystalWidth/2,
                 crystalWidth/2,
                 crystalLength/2);   
                           // size
  fCrystalLog 
    = new G4LogicalVolume(crystal_box,    // its solid
                          pbWO4,          // its material
                          "CrystalLog"); // its name
          
  for (G4int i=0; i<nbOfCrystals; i++)
    {
      
      for (G4int j=0; j<nbOfCrystals; j++)
        {  
          G4int n =  i*10+j;
          G4ThreeVector crystalPos((i*crystalWidth)-135,
                                   (j*crystalWidth)-135,0 );
          fCrystalPhys[n] 
            = new G4PVPlacement(0,               // no rotation
                                crystalPos,      // translation
                                fCrystalLog,
                                "crystal",       // its name
                                caloLog,
                                false,
                                i);
        }
    }  
  G4cout << "There are " << nbOfCrystals <<
    " crystals per row in the calorimeter, so in total "<<
    nbOfCrystals*nbOfCrystals << " crystals" << G4endl;  
  G4cout << "They have width of  " << crystalWidth /cm <<
    "  cm and a length of  " <<  crystalLength /cm
         <<" cm. The Material is "<< pbWO4 << G4endl;
  
  
  experimentalHallLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* caloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  caloLog->SetVisAttributes(caloVisAtt);
  fCrystalLog->SetVisAttributes(crystalVisAtt);


  // Parametrisation region defined in ExGflash2ParallelWorld
  
  return experimentalHallPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflash2DetectorConstruction::ConstructSDandField()
{
  // -- sensitive detectors:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  ExGflash2SensitiveDetector* CaloSD
    = new ExGflash2SensitiveDetector("Calorimeter",this);
  SDman->AddNewDetector(CaloSD);
  fCrystalLog->SetSensitiveDetector(CaloSD);
  
  // Fast simulation implemented in ExGflash2ParallelWorld
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
