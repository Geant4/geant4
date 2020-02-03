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
/// \file persistency/P03/src/ExTGDetectorConstruction.cc
/// \brief Implementation of the ExTGDetectorConstruction class

#include "G4tgbVolumeMgr.hh"
#include "ExTGDetectorConstruction.hh"
#include "G4tgrMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstruction::ExTGDetectorConstruction()
{ 
  fMessenger = new G4tgrMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstruction::~ExTGDetectorConstruction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExTGDetectorConstruction::Construct()
{
  //------------------------------------------------ 
  // Define one or several text files containing the geometry description
  //------------------------------------------------ 
  G4String filename = "g4geom.txt";
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);

  //------------------------------------------------ 
  // Read the text files and construct the GEANT4 geometry
  //------------------------------------------------ 
  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();

  return physiWorld;
}
