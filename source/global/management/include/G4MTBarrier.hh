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
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// Class Description:
//
// This class defines a synchronization point between threads: a master
// and a pool of workers.
// A barrier is a (shared) instance of this class. Master sets the number
// of active threads to wait for, then it waits for workers to become ready
// calling the method WaitForReadyWorkers(). The master thread will block on this
// call.
// Each of the workers calls ThisWorkerReady() when it is ready to continue.
// It will block on this call.
// When all worker threads have called ThisWorkerReady and are waiting the
// master will release the barrier and execution will continue.
//
// User code can implement more advanced barriers that require exchange
// of a message between master and threads inheriting from this class as in:
//    class Derived : public G4MTBarrier {
//       G4Mutex mutexForMessage;
//       SomeType message;
//       void MethodCalledByWorkers() {
//            G4MTBarrirer::ThisWorkerReady();
//            G4AutoLock l(&mutexForMessage);
//            [... process message ...]
//       }
//      void WaitForReadyWorkers() override {
//             Wait(); <== Mandatory
//             [.. process message ...] <== User code between the two calls
//             ReleaseBarrier(); <== Mandatory
//      }
//      void MethodCalledByMaster() { WaitForReadyWorkers(); }
//    }
// User code can also achieve the same results as before using the granular
// methods LoopWaitingWorkers and ResetCounterAndBroadcast methods in the
// master. For examples of usage of this class see G4MTRunManager
//
// G4MTBarrier.hh
//
//  Created on: Feb 10, 2016
//      Author: adotti
//
// =====================================
// Barriers mechanism
// =====================================
// We want to implement barriers.
// We define a barrier has a point in which threads synchronize.
// When workers threads reach a barrier they wait for the master thread a
// signal that they can continue. The master thread broadcast this signal
// only when all worker threads have reached this point.
// Currently only three points require this sync in the life-time of a G4 applicattion:
// Just before and just after the for-loop controlling the thread event-loop.
// Between runs.
//
// The basic algorithm of each barrier works like this:
// In the master:
//   WaitWorkers()  {
//    while (true)
//    {
//     G4AutoLock l(&counterMutex);                          || Mutex is locked (1)
//     if ( counter == nActiveThreads ) break;
//     G4CONDITIONWAIT( &conditionOnCounter, &counterMutex); || Mutex is atomically released and wait, upon return locked (2)
//    }                                                      || unlock mutex
//    G4AutoLock l(&counterMutex);                           || lock again mutex (3)
//    G4CONDITIONBROADCAST( &doSomethingCanStart );          || Here mutex is locked (4)
//   }                                                       || final unlock (5)
// In the workers:
//   WaitSignalFromMaster() {
//    G4AutoLock l(&counterMutex);                           || (6)
//    ++counter;
//    G4CONDITIONBROADCAST(&conditionOnCounter);             || (7)
//    G4CONDITIONWAIT( &doSomethingCanStart , &counterMutex);|| (8)
//   }
// Each barriers requires 2 conditions and one mutex, plus a counter.
// Important note: the thread calling broadcast should hold the mutex
// before calling broadcast to obtain predictible behavior
// http://pubs.opengroup.org/onlinepubs/7908799/xsh/pthread_cond_broadcast.html
// Also remember that the wait for condition will atomically release the mutex
// and wait on condition, but it will lock again on mutex when returning
// Here it is how the control flows.
// Imagine master starts and only one worker (nActiveThreads==1)
// Master       |    Worker        | counter | Who holds mutex
// Gets to (1)  |   Blocks on (6)  | 0       | M
// Waits in (2) |                  | 0       | -
//              |  Arrives to (7)  | 1       | W
//              |  Waits in (8)    | 1       | -
// Gets to (1)  |                  | 1       | M
// Jumps to (3) |                  | 1       | M
// End          |                  | 1       | -
//              | End              | 1       | -
// Similarly for more than one worker threads or if worker starts

#ifndef G4MTBARRIER_HH_
#define G4MTBARRIER_HH_
#include "G4Threading.hh"

#ifdef WIN32
#include "windefs.hh"
#endif

class G4MTBarrier
{
public:
  G4MTBarrier() : G4MTBarrier(1) {}
  virtual ~G4MTBarrier() {}
  G4MTBarrier(const G4MTBarrier&) = delete;
  G4MTBarrier& operator=(const G4MTBarrier&) = delete;
  //on explicitly defaulted move at
  //https://msdn.microsoft.com/en-us/library/dn457344.aspx
  //G4MTBarrier(G4MTBarrier&&) = default;
  //G4MTBarrier& operator=(G4MTBarrier&&)  = default;
  G4MTBarrier( unsigned int numThreads );
  void ThisWorkerReady();
  virtual void WaitForReadyWorkers();
  inline void SetActiveThreads( unsigned int val ) { m_numActiveThreads = val; }
  void ResetCounter();
  unsigned int GetCounter();
  void Wait();
  void ReleaseBarrier();
  inline void Wait( unsigned int numt ) {
    SetActiveThreads( numt );
    Wait();
  }
private:
  unsigned int m_numActiveThreads;
  unsigned int m_counter;
  G4Mutex m_mutex;
  G4Condition m_counterChanged;
  G4Condition m_continue;
#if defined(WIN32)
  CRITICAL_SECTION cs1;
  CRITICAL_SECTION cs2;
#endif

};

#endif /* G4MTBARRIER_HH_ */
