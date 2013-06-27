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
#include "G4VtbbJob.hh"

#include "G4AutoLock.hh"
#include "G4UImanager.hh"
#include "G4tbbRunManager.hh"

//Global Mutex
G4Mutex init_worker_mutex = G4MUTEX_INITIALIZER;

G4VtbbJob::G4VtbbJob(const G4String& macro) : macroFile(macro)
{
}

G4VtbbJob::~G4VtbbJob() {
}

void G4VtbbJob::ThreadSafeInitSetup(G4tbbRunManager* rm)
{
  //Acquire global lock to serialize initialization of worker thread
  G4AutoLock lock(&init_worker_mutex);
  //Step 1-Now instead of calling UserActions, call the worker-thread equivalent
  this->InitSetup(rm);
  //Step 2- Now Do whatever is needed to finalize the job
  this->JobPrepare(rm);
  //Step 3- Process macro file
  //If Macro File contains a /run/beamOn command simulation should not
  //be used: use directly the G4RunManager::BeamOn(...) cmd
  rm->Initialize();

  if ( macroFile != "" ) {
    G4UImanager * UI = G4UImanager::GetUIpointer();  
    G4String command = "/control/execute ";
    UI->ApplyCommand(command+macroFile);
  }
}

void G4VtbbJob::InitRun( G4tbbRunManager* rm)
{
  //This is called once by the master thread to execute the macro file
  //Step 0- The G4tbbRunManager instance has been created by the caller
  if ( rm == NULL ) {
    //Throw exception!
    G4Exception("G4VtbbJob::InitRun()", "RunTbbJob001",
                FatalException, "G4RunManager is not defined!");
  }
  //Step 1- Now call User actions definitions and creation

  //  this->CreateDetector(rm); - Replace the line below: better name
  this->UserActions(rm);
  //Step 2- Now Do whatever is needed to finalize the job
  this->JobPrepare(rm);

  //Step 3- Start simulation
  //If Macro File contains a /run/beamOn command simulation should not
  //be used: use directly the G4RunManager::BeamOn(...) cmd
  rm->Initialize();

  if ( macroFile != "" ) {
    G4UImanager * UI = G4UImanager::GetUIpointer();  
    G4String command = "/control/execute ";
    UI->ApplyCommand(command+macroFile);
  }
}
