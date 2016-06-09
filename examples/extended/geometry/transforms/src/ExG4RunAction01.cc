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
// $Id$
//
/// \file ExG4RunAction01.cc
/// \brief Implementation of the ExG4RunAction01 class

#include "ExG4RunAction01.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
//#include "G4ios.hh"
//#include <iomanip>

#include "Randomize.hh"

const G4String ExG4RunAction01::fgkDefaultRndmFileName = "run0.rndm";

//_____________________________________________________________________________
ExG4RunAction01::ExG4RunAction01()
  : G4UserRunAction(),
    fMessenger(this),
    fRndmFileName(fgkDefaultRndmFileName),
    fVerboseLevel(1),
    fSaveRndm(false),
    fReadRndm(false),
    fAutoSeed(false)
{
/// Standard constructor
}

//_____________________________________________________________________________
ExG4RunAction01::~ExG4RunAction01()
{
/// Destructor
}

//_____________________________________________________________________________
void ExG4RunAction01::BeginOfRunAction(const G4Run* run)
{ 
  if ( fVerboseLevel > 0 ) { 
    G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  }  
  
  // read Rndm status
  if ( fReadRndm ) {
    G4cout << "\n---> rndm status restored from file: " << fRndmFileName << G4endl;
    CLHEP::HepRandom::restoreEngineStatus(fRndmFileName);
    CLHEP::HepRandom::showEngineStatus();
  }   

  // automatic (time-based) random seeds for each run
  if ( fAutoSeed ) {
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  
    G4cout << "*******************" << G4endl;
    G4cout << "*** AUTOSEED ON ***" << G4endl;
    G4cout << "*******************" << G4endl;
    
    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    CLHEP::HepRandom::setTheSeeds(seeds);
    CLHEP::HepRandom::showEngineStatus();
  }
  
  // save Rndm status
  if ( fSaveRndm ) { 
    G4int runNumber = run->GetRunID();
    std::ostringstream fileName;
    fileName << "run" << runNumber << ".rndm";
    CLHEP::HepRandom::saveEngineStatus(fileName.str().c_str()); 
    CLHEP::HepRandom::showEngineStatus();
  }  
  
  if ( G4VVisManager::GetConcreteInstance() )  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
  }  
}

//_____________________________________________________________________________
void ExG4RunAction01::EndOfRunAction(const G4Run* run)
{
  // save Rndm status
  if ( fSaveRndm ) { 
    G4int runNumber = run->GetRunID();
    std::ostringstream fileName;
    fileName << "run" << runNumber << "end.rndm";
    CLHEP::HepRandom::saveEngineStatus(fileName.str().c_str()); 
    CLHEP::HepRandom::showEngineStatus();
  }     

  if ( G4VVisManager::GetConcreteInstance() ) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
}
