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
// Author: Jonathan Madsen (May 28st 2020)
//
// class description:
//
//      This class is used for multi-threaded Geant4.
//      It encapsulates the mechanism of starting/stopping threads.

#ifndef G4UserTaskThreadInitialization_hh
#define G4UserTaskThreadInitialization_hh

class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4WorkerThread;
class G4WorkerTaskRunManager;

#include "G4Threading.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "Randomize.hh"

class G4UserTaskThreadInitialization : public G4UserWorkerThreadInitialization
{
 public:  // with description
  G4UserTaskThreadInitialization();
  virtual ~G4UserTaskThreadInitialization();

  virtual G4Thread* CreateAndStartWorker(G4WorkerThread* workerThreadContext);
  //  Called by the kernel to create a new thread/worker
  //  and start work.
  // Usere should not re-implement this function (in derived class), except only
  // if he/she wants to verwrite the default threading model (see StartThread
  // function)

  virtual void SetupRNGEngine(const CLHEP::HepRandomEngine* aRNGEngine) const;
  // Called by worker threads to set the Random Number Generator Engine
  // The default implementation "clones" the engine from the master thread
  // User needs to re-implement this method if using a non-standard
  // RNG Engine (i.e. a different one w.r.t. the one provided in the CLHEP
  // version supported by G4.
  // Important: this method is called by all threads at the same time
  //   if is user responsibilitiy to make it thread-safe

  virtual void JoinWorker(G4Thread* aThread);
  // Called by the kernel when threads need to be terminated. Implements logic
  // of joining the aThread. Calling thread will wait for aThread to end. Usere
  // should not re-implement this function (in derived class), except only if
  // he/she wants to verwrite the default threading model (see StartThread
  // function)

  virtual G4WorkerRunManager* CreateWorkerRunManager() const;
  // Called by StartThread function to create a run-manager implementing worker
  // behvior. User should re-implemtn this function in derived class to
  // instantiate his/her user-defined WorkerRunManager. By default this method
  // instantiates G4WorkerTaskRunManager object.
};

#endif  // G4UserTaskThreadInitialization_hh
