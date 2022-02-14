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
// G4UserWorkerThreadInitialization implementation
//
// Authors: M.Asai, A.Dotti (SLAC), 16 September 2013
// --------------------------------------------------------------------

#include <sstream>

#include "G4UserWorkerThreadInitialization.hh"
#include "G4AutoLock.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4UImanager.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4WorkerRunManager.hh"
#include "G4WorkerThread.hh"

// Will need this for TPMalloc
//#ifdef G4MULTITHREADED
//#define TPMALLOCDEFINESTUB
//#include "tpmalloc/tpmallocstub.h"
//#endif

// --------------------------------------------------------------------
#ifdef G4MULTITHREADED
G4Thread* G4UserWorkerThreadInitialization::
CreateAndStartWorker(G4WorkerThread* wTC)
{
  // Note: this method is called by G4MTRunManager,
  // here we are still sequential.
  // Create a new thread/worker structure
  G4Thread* worker = new G4Thread;
  G4THREADCREATE(worker, &G4MTRunManagerKernel::StartThread, wTC);
  return worker;
}
#else
G4Thread* G4UserWorkerThreadInitialization::
CreateAndStartWorker(G4WorkerThread*)
{
  return new G4Thread;
}
#endif

// --------------------------------------------------------------------
#ifdef G4MULTITHREADED
void G4UserWorkerThreadInitialization::JoinWorker(G4Thread* aThread)
{
  G4THREADJOIN(*aThread);
}
#else  // Avoid compilation warning in sequential
void G4UserWorkerThreadInitialization::JoinWorker(G4Thread*)
{
}
#endif

// --------------------------------------------------------------------
G4UserWorkerThreadInitialization::G4UserWorkerThreadInitialization()
{
}

// --------------------------------------------------------------------
G4UserWorkerThreadInitialization::~G4UserWorkerThreadInitialization()
{
}

// --------------------------------------------------------------------
namespace
{
  G4Mutex rngCreateMutex = G4MUTEX_INITIALIZER;
}

// --------------------------------------------------------------------
void G4UserWorkerThreadInitialization::
SetupRNGEngine(const CLHEP::HepRandomEngine* aNewRNG) const
{
  G4AutoLock l(&rngCreateMutex);
  // No default available, let's create the instance of random stuff
  // A Call to this just forces the creation to defaults
  G4Random::getTheEngine();
  // Poor man's solution to check which RNG Engine is used in master thread
  CLHEP::HepRandomEngine* retRNG = nullptr;

  // Need to make these calls thread safe
  if(dynamic_cast<const CLHEP::HepJamesRandom*>(aNewRNG))
  {
    retRNG = new CLHEP::HepJamesRandom;
  }
  if(dynamic_cast<const CLHEP::MixMaxRng*>(aNewRNG))
  {
    retRNG = new CLHEP::MixMaxRng;
  }
  if(dynamic_cast<const CLHEP::RanecuEngine*>(aNewRNG))
  {
    retRNG = new CLHEP::RanecuEngine;
  }
  if(dynamic_cast<const CLHEP::RanluxppEngine*>(aNewRNG))
  {
    retRNG = new CLHEP::RanluxppEngine;
  }
  if(dynamic_cast<const CLHEP::Ranlux64Engine*>(aNewRNG))
  {
    const CLHEP::Ranlux64Engine* theRNG =
      dynamic_cast<const CLHEP::Ranlux64Engine*>(aNewRNG);
    retRNG = new CLHEP::Ranlux64Engine(123, theRNG->getLuxury());
  }
  if(dynamic_cast<const CLHEP::MTwistEngine*>(aNewRNG))
  {
    retRNG = new CLHEP::MTwistEngine;
  }
  if(dynamic_cast<const CLHEP::DualRand*>(aNewRNG))
  {
    retRNG = new CLHEP::DualRand;
  }
  if(dynamic_cast<const CLHEP::RanluxEngine*>(aNewRNG))
  {
    const CLHEP::RanluxEngine* theRNG =
      dynamic_cast<const CLHEP::RanluxEngine*>(aNewRNG);
    retRNG = new CLHEP::RanluxEngine(123, theRNG->getLuxury());
  }
  if(dynamic_cast<const CLHEP::RanshiEngine*>(aNewRNG))
  {
    retRNG = new CLHEP::RanshiEngine;
  }

  if(retRNG != nullptr)
  {
    G4Random::setTheEngine(retRNG);
  }
  else
  {
    // Does a new method, such as aNewRng->newEngine() exist to clone it ?
    G4ExceptionDescription msg;
    msg << " Unknown type of RNG Engine - " << G4endl
        << " Can cope only with HepJamesRandom, MixMaxRng, Ranecu, Ranlux64,"
        << " Ranlux++, MTwistEngine, DualRand, Ranlux or Ranshi." << G4endl
        << " Cannot clone this type of RNG engine, as required for this thread"
        << G4endl << " Aborting " << G4endl;
    G4Exception("G4UserWorkerThreadInitialization::SetupRNGEngine()",
                "Run0122", FatalException, msg);
  }
}

// --------------------------------------------------------------------
G4WorkerRunManager*
G4UserWorkerThreadInitialization::CreateWorkerRunManager() const
{
  return new G4WorkerRunManager();
}
