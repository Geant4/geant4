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
/// \file ExG4EventAction01.cc
/// \brief Implementation of the ExG4EventAction01 class

#include "ExG4EventAction01.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"

#include <sstream>

const  G4int ExG4EventAction01::fgkDefaultVerboseLevel = 1;
const  G4int ExG4EventAction01::fgkDefaultPrintModulo = 10000;

//_____________________________________________________________________________
ExG4EventAction01::ExG4EventAction01()
 : G4UserEventAction(),
   fMessenger(this),
   fVerboseLevel(fgkDefaultVerboseLevel),
   fPrintModulo(fgkDefaultPrintModulo),
   fSaveRndm(false)
{
/// Standard constructor
}

//_____________________________________________________________________________
ExG4EventAction01::~ExG4EventAction01()
{
/// Destructor
}

//_____________________________________________________________________________
void ExG4EventAction01::BeginOfEventAction(const G4Event* event)
{
  // Print event info
  //
  G4int eventNumber = event->GetEventID();
  if ( eventNumber%fPrintModulo == 0 ) {
    G4cout << "\n---> Begin of Event: " << eventNumber << G4endl;
  }   
     
  // Print verbose info
  //
  if ( fVerboseLevel > 1 ) {
    G4cout << "<<< Event  " << eventNumber << " started." << G4endl;
  }   
}

//_____________________________________________________________________________
void ExG4EventAction01::EndOfEventAction(const G4Event* event)
{
  // Print verbose info
  //
  if ( fVerboseLevel > 1 ) {
    G4cout << "<<< Event  " << event->GetEventID() << " ended." << G4endl;
  }  
  
  // Save rndm status
  //
  if ( fSaveRndm ) {     
    const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
    G4int runNumber = run->GetRunID();
    G4int eventNumber = event->GetEventID();
    std::ostringstream fileName;
    fileName << "run" << runNumber << "event" << eventNumber << ".rndm";
    CLHEP::HepRandom::saveEngineStatus(fileName.str().c_str()); 

    if ( eventNumber%fPrintModulo == 0 ) {
      G4cout << "\n---> End of Event: " << eventNumber << G4endl;
      CLHEP::HepRandom::showEngineStatus();
    }
  }       
}
