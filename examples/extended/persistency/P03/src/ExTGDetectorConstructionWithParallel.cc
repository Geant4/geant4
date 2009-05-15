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
// $Id: ExTGDetectorConstructionWithParallel.cc,v 1.1 2009-05-15 16:39:14 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------------

#include "G4tgbVolumeMgr.hh"
#include "ExTGDetectorConstructionWithParallel.hh"
#include "G4tgrMessenger.hh"
#include "G4tgbParallelGeomMgr.hh"
#include "G4tgbDetectorBuilder.hh"


#include "G4VUserParallelWorld.hh"

// ---------------------------------------------------------------------------
ExTGDetectorConstructionWithParallel::ExTGDetectorConstructionWithParallel()
{
  messenger = new G4tgrMessenger;
}

// ---------------------------------------------------------------------------
ExTGDetectorConstructionWithParallel::~ExTGDetectorConstructionWithParallel()
{
  delete messenger;
}

// ---------------------------------------------------------------------------
void ExTGDetectorConstructionWithParallel::ConstructParallelWorlds()
{
  //------------------------------------------------ 
  // Define one or several text files containing the geometry description
  //------------------------------------------------ 
  G4String filename = "g4geom.txt";
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);

  G4String filenameP = "g4geom_parallel.txt";
  volmgr->AddTextFileParallel(filenameP,1);

  //------------------------------------------------ 
  // Read the text files 
  //------------------------------------------------ 
  //-  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();
  G4tgbDetectorBuilder* detectorBuilder = volmgr->GetDetectorBuilder();
  theTgrVoltop = detectorBuilder->ReadDetector();

  //------------------------------------------------ 
  // Construct the Geant4 parallel worlds
  //------------------------------------------------ 
  std::vector<G4VUserParallelWorld*> parallelWorlds = G4tgbParallelGeomMgr::GetInstance()->CreateParalleWorlds();
  for( size_t ii = 0; ii < parallelWorlds.size(); ii++ ) {
    G4cout << " RegisterParallelWorld " << parallelWorlds[ii]->GetName() << G4endl;
    RegisterParallelWorld( parallelWorlds[ii] );
   
  }
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* ExTGDetectorConstructionWithParallel::Construct()
{
  ConstructParallelWorlds();

  //------------------------------------------------ 
  // Construct the Geant4 mass geometry
  //------------------------------------------------ 
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  G4tgbDetectorBuilder* detectorBuilder = volmgr->GetDetectorBuilder();
  G4VPhysicalVolume* physiWorld = detectorBuilder->ConstructDetector(theTgrVoltop);

  return physiWorld;
}
