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
//
/// \file persistency/P03/src/ExTGDetectorConstructionWithSD.cc
/// \brief Implementation of the ExTGDetectorConstructionWithSD class

#include "ExTGDetectorConstructionWithSD.hh"
#include "ExTGTrackerSD.hh"

#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

#include "G4tgbVolumeMgr.hh"
#include "G4tgrMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstructionWithSD::ExTGDetectorConstructionWithSD()
{
  fMessenger = new G4tgrMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstructionWithSD::~ExTGDetectorConstructionWithSD()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExTGDetectorConstructionWithSD::Construct()
{
  //------------------------------------------------ 
  // Define one or several text files containing the geometry description
  //------------------------------------------------ 
  G4String filename = "g4geom_SD.txt";
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);

  //------------------------------------------------ 
  // Read the text files and construct the GEANT4 geometry
  //------------------------------------------------ 
  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();

  return physiWorld;
}

void ExTGDetectorConstructionWithSD::ConstructSDandField()
{
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4String trackerChamberSDname = "ExTextGeom/TrackerChamberSD";
  ExTGTrackerSD* aTrackerSD = new ExTGTrackerSD( trackerChamberSDname );
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  G4LogicalVolume * logicChamber =
    G4tgbVolumeMgr::GetInstance()->FindG4LogVol("Chamber",0);
  if(logicChamber)
  {
    SetSensitiveDetector("Chamber", aTrackerSD );
  }
  else
  {
    G4Exception("ExTGDetectorConstructionWithSD::Construct()",
                "InvalidGeometry", JustWarning,
                "Volume does not exists in geometry: Chamber");
  }

}
