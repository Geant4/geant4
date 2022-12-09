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
/// \file ExGflash3ParallelWorld.cc
/// \brief Implementation of the ExGflash3ParallelWorld class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// User Classes
#include "ExGflash3ParallelWorld.hh"
#include "ExGflash3SensitiveDetector.hh"

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

ExGflash3ParallelWorld::ExGflash3ParallelWorld(G4String aWorldName)
  :G4VUserParallelWorld(aWorldName), fCrystalLog(nullptr), fCrystalPhys{}
{
  G4cout<<"ExGflash3ParallelWorld::Parralel world constructor"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflash3ParallelWorld::~ExGflash3ParallelWorld() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflash3ParallelWorld::Construct()
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
  G4LogicalVolume* caloLog 
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
                    caloLog,
                    "calorimeter",
                    ghostLogicalVolume,
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
                          dummy,          // its material
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
         <<" cm. " << G4endl;
  

  G4VisAttributes* caloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  caloLog->SetVisAttributes(caloVisAtt);
  fCrystalLog->SetVisAttributes(crystalVisAtt);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflash3ParallelWorld::ConstructSD()
{
  // -- sensitive detectors:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  ExGflash3SensitiveDetector* CaloSD
    = new ExGflash3SensitiveDetector("Calorimeter",this);
  SDman->AddNewDetector(CaloSD);
  fCrystalLog->SetSensitiveDetector(CaloSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
