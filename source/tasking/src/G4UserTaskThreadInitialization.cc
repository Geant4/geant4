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

#include "G4UserTaskThreadInitialization.hh"
#include "G4AutoLock.hh"
#include "G4TaskRunManagerKernel.hh"
#include "G4UImanager.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4WorkerRunManager.hh"
#include "G4WorkerTaskRunManager.hh"
#include "G4WorkerThread.hh"
#include "globals.hh"
#include <sstream>

//============================================================================//

namespace
{
  G4Mutex rngCreateMutex;
}

//============================================================================//

G4Thread* G4UserTaskThreadInitialization::CreateAndStartWorker(G4WorkerThread*)
{
  // Note: this method is called by G4MTRunManager, here we are still sequential
  // Create a new thread/worker structure
  return nullptr;
}

//============================================================================//

// Avoid compilation warning in sequential
void G4UserTaskThreadInitialization::JoinWorker(G4Thread* aThread)
{
  if(aThread)
  {
    G4THREADJOIN(*aThread);
  }
}

//============================================================================//

G4UserTaskThreadInitialization::G4UserTaskThreadInitialization() {}

//============================================================================//

G4UserTaskThreadInitialization::~G4UserTaskThreadInitialization() {}

//============================================================================//

void G4UserTaskThreadInitialization::SetupRNGEngine(
  const CLHEP::HepRandomEngine* aNewRNG) const
{
  G4AutoLock l(&rngCreateMutex);
  // No default available, let's create the instance of random stuff
  // A Call to this just forces the creation to defaults
  G4Random::getTheEngine();
  // Poor man's solution to check which RNG Engine is used in master thread
  CLHEP::HepRandomEngine* retRNG = nullptr;

  // Need to make these calls thread safe
  if(dynamic_cast<const CLHEP::HepJamesRandom*>(aNewRNG))
    retRNG = new CLHEP::HepJamesRandom;
  if(dynamic_cast<const CLHEP::MixMaxRng*>(aNewRNG))
    retRNG = new CLHEP::MixMaxRng;
  if(dynamic_cast<const CLHEP::RanecuEngine*>(aNewRNG))
    retRNG = new CLHEP::RanecuEngine;
  if(dynamic_cast<const CLHEP::Ranlux64Engine*>(aNewRNG))
    retRNG = new CLHEP::Ranlux64Engine;
  if(dynamic_cast<const CLHEP::RanluxppEngine*>(aNewRNG))
    retRNG = new CLHEP::RanluxppEngine;
  if(dynamic_cast<const CLHEP::MTwistEngine*>(aNewRNG))
    retRNG = new CLHEP::MTwistEngine;
  if(dynamic_cast<const CLHEP::DualRand*>(aNewRNG))
    retRNG = new CLHEP::DualRand;
  if(dynamic_cast<const CLHEP::RanluxEngine*>(aNewRNG))
    retRNG = new CLHEP::RanluxEngine;
  if(dynamic_cast<const CLHEP::RanshiEngine*>(aNewRNG))
    retRNG = new CLHEP::RanshiEngine;

  if(retRNG != nullptr)
    G4Random::setTheEngine(retRNG);
  else
  {
    // Does a new method, such as aNewRng->newEngine() exist to clone it ?
    G4ExceptionDescription msg;
    msg << " Unknown type of RNG Engine - " << G4endl
        << " Can cope only with HepJamesRandom, MixMaxRng, Ranecu, Ranlux64,"
        << " Ranlux++, MTwistEngine, DualRand, Ranlux or Ranshi." << G4endl
        << " Cannot clone this type of RNG engine, as required for this thread"
        << G4endl << " Aborting... " << G4endl;
    G4Exception("G4UserTaskInitializition::SetupRNGEngine()", "Run0122",
                FatalException, msg);
  }
}

//============================================================================//

G4WorkerRunManager* G4UserTaskThreadInitialization::CreateWorkerRunManager()
  const
{
  return new G4WorkerTaskRunManager();
}

//============================================================================//
