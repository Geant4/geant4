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
/// \file persistency/P03/src/ExTGDetectorConstructionWithCuts.cc
/// \brief Implementation of the ExTGDetectorConstructionWithCuts class

#include "ExTGDetectorConstructionWithCuts.hh"
#include "ExTGRCDetectorBuilder.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgrMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstructionWithCuts::ExTGDetectorConstructionWithCuts()
{
  fMessenger = new G4tgrMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExTGDetectorConstructionWithCuts::~ExTGDetectorConstructionWithCuts()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExTGDetectorConstructionWithCuts::Construct()
{
  //------------------------------------------------ 
  // Define one or several text files containing the geometry description
  //------------------------------------------------ 
  G4String filename = "g4geom_cutsPerRegion.txt";
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);

  //------------------------------------------------ 
  // Use your own detector builder, that will invoke your own line processor
  //------------------------------------------------ 
  ExTGRCDetectorBuilder* gtb = new ExTGRCDetectorBuilder;
  volmgr->SetDetectorBuilder( gtb );

  const G4tgrVolume* tgrVoltop = gtb->ReadDetector();
  G4VPhysicalVolume* physiWorld = gtb->ConstructDetector(tgrVoltop);

  //  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();

  return physiWorld;
}
