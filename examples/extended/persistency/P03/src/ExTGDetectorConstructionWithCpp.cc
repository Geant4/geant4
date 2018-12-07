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
/// \file persistency/P03/src/ExTGDetectorConstructionWithCpp.cc
/// \brief Implementation of the ExTGDetectorConstructionWithCpp class

#include "G4tgbVolumeMgr.hh"
#include "ExTGDetectorConstructionWithCpp.hh"
#include "G4tgrMessenger.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstructionWithCpp::ExTGDetectorConstructionWithCpp()
{
  fMessenger = new G4tgrMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstructionWithCpp::~ExTGDetectorConstructionWithCpp()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExTGDetectorConstructionWithCpp::Construct()
{
  //------------------------------------------------ 
  // Define one or several text files containing the geometry description
  //------------------------------------------------ 
  G4String filename = "g4geom_simple.txt";
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);

  //------------------------------------------------ 
  // Read the text files and construct the GEANT4 geometry
  //------------------------------------------------ 
  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();

  //------------------------------------------------ 
  // Build another volume and place it in the text geometry
  //------------------------------------------------ 
  G4Box* solid = new G4Box("mybox",2.,2.,3.);
  
  G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");

  G4LogicalVolume* logicVol = new G4LogicalVolume(solid,water,"mybox",0,0,0);

  //----- Place the volume in the world of the text geometry
  G4LogicalVolume* textWorldVol = physiWorld->GetLogicalVolume();

  // NOTE: if you want to place your full text geometry in the C++ geometry,
  //       you should use this G4LogicalVolume and place it in a C++ logical
  //       volume
  new G4PVPlacement(0,               // no rotation
                    G4ThreeVector(0.,-20.,0),  // at (x,y,z)
                    logicVol,        // its logical volume    
                    "myBox1",        // its name
                    textWorldVol,    // its mother  volume
                    false,           // no boolean operations
                    0);              // copy number 
  
  //----- Place the volume inside a volume of the text geometry
  G4LogicalVolume* textSphereVol =
    G4tgbVolumeMgr::GetInstance()->FindG4LogVol("sphere",1);
    new G4PVPlacement(0,             // no rotation
                    G4ThreeVector(0.,0.,0),  // at (x,y,z)
                    logicVol,        // its logical volume    
                    "myBox2",        // its name
                    textSphereVol,   // its mother  volume
                    false,           // no boolean operations
                    1);              // copy number 
  
  return physiWorld;
}
