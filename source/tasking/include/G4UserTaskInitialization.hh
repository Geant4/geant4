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
//  Author: Jonathan Madsen (May 28st 2020)
//
// class description:
//
//      This class is used for multi-threaded Geant4.
//      The object of this class can be set to G4MTRunManager, but not to
//      G4RunManager. G4UserTaskInitialization class has five virtual methods as
//      the user hooks which are invoked at several occasions of the life cycle
//      of each thread.
//
//      - virtual void WorkerInitialize() const
//       This method is called after the tread is created but before the
//       G4WorkerTaskRunManager is instantiated.
//      - virtual void WorkerStart() const
//       This method is called once at the beginning of simulation job when
//       kernel classes and user action classes have already instantiated but
//       geometry and physics have not been yet initialized. This situation is
//       identical to "PreInit" state in the sequential mode.
//      - virtual void WorkerRunStart() const
//       This method is called before an event loop. Geometry and physics have
//       already been set up for the thread. All threads are synchronized and
//       ready to start the local event loop. This situation is identical to
//       "Idle" state in the sequential mode.
//      - virtual void WorkerRunEnd() const
//       This method is called for each thread when the local event loop is
//       done, but before the synchronization over threads.
//      - virtual void WorkerStop() const
//       This method is called once at the end of simulation job.
//
//      Note: This object should be instantiated only once and set to
//      G4MTRunManager,
//       while these five methods are invoked for each worker thread. Thus, to
//       store thread-local objects, use G4ThreadLocal keyword.
//

#ifndef G4UserTaskInitialization_hh
#define G4UserTaskInitialization_hh

#include "G4UserWorkerInitialization.hh"

class G4UserTaskInitialization : public G4UserWorkerInitialization
{
 public:  // with description
  G4UserTaskInitialization();
  virtual ~G4UserTaskInitialization();

  virtual void WorkerInitialize() const;
  // This method is called after the tread is created but before the
  // G4WorkerTaskRunManager is instantiated.

  virtual void WorkerStart() const;
  // This method is called once at the beginning of simulation job
  // when kernel classes and user action classes have already instantiated
  // but geometry and physics have not been yet initialized. This situation
  // is identical to "PreInit" state in the sequential mode.

  virtual void WorkerRunStart() const;
  // This method is called before an event loop. Geometry and physics have
  // already been set up for the thread. All threads are synchronized and
  // ready to start the local event loop. This situation is identical to
  // "Idle" state in the sequential mode.

  virtual void WorkerRunEnd() const;
  // This method is called for each thread, when the local event loop has
  // finished but before the synchronization over threads.

  virtual void WorkerStop() const;
  // This method is called once at the end of simulation job.
  // Implement here a clean up action.
};

#endif  // G4UserTaskInitialization_hh
