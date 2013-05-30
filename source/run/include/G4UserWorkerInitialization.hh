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
//  25 Feb 2013: Andrea Dotti, first implementation
//
// class description:
//
//      This class is used for multi-threaded Geant4
//      It encapsulates the user-defined instantiations of objects
//      and setup preparation for each worker-thread.
//      In its most simple implementation, the virtual method:
//      WorkerStart() implements what in a sequential G4 application
//      is done in the "main" function.
//      For example a MT application can be defined with:
//          class myWorkerInitialization : public G4UserWorkerInitialization {
//          public:
//              void WorkerStart() {
//                  SetUserAction( //…//);
//                  //other usercode
//              }
//          };
//
//          int main(int, char**) {
//              G4MTRunManager* rm = new G4MTRunManager;
//              rm->SetUserInitialization( new myWorkerInitialization );
//              rm->SetNumberThreads( nThreads );
//              rm->BeamOn( … );
//          }

#ifndef G4UserWorkerInitialization_hh
#define G4UserWorkerInitialization_hh

class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

#include "G4Threading.hh"
#include "G4WorkerThread.hh"
#include "G4WorkerRunManager.hh"
#include "G4MTRunManager.hh"

class G4UserWorkerInitialization {
public: // with description
    G4UserWorkerInitialization();
    virtual ~G4UserWorkerInitialization();

    virtual void WorkerInitialize() const;
    // This method is called after the tread is created but before the
    // G4WorkerRunManager is instantiated.

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
 
    virtual G4Thread* CreateAndStartWorker(G4WorkerThread* workerThreadContext);
    //  Called by the kernel to create a new thread/worker
    //  and start work.
    // Usere should not re-implement this function (in derived class), except only if he/she
    // wants to verwrite the default threading model (see StartThread function)

    virtual void SetupRNGEngine(const CLHEP::HepRandomEngine* aRNGEngine) const;
    // Called by worker threads to set the Random Number Generator Engine
    // The default implementation "clones" the engine from the master thread
    // User needs to re-implement this method if using a non-standard
    // RNG Engine (i.e. a different one w.r.t. the one provided in the CLHEP
    // version supported by G4.
    // Important: this method is called by all threads at the same time
    //   if is user responsibilitiy to make it thread-safe

protected:   
    static G4ThreadLocal G4WorkerThread* wThreadContext;
    
protected: // with description
    void SetUserAction(G4VUserPrimaryGeneratorAction*) const;
    void SetUserAction(G4UserRunAction*) const;
    void SetUserAction(G4UserEventAction*) const;
    void SetUserAction(G4UserStackingAction*) const;
    void SetUserAction(G4UserTrackingAction*) const;
    void SetUserAction(G4UserSteppingAction*) const;
    // These methods should be used to define user's action classes.

private:
    static void* StartThread( void* context);
    // This static function is used to start a pthread-based
    // worker.
    // context is the instance of type G4UserWorkerAction.
    // Method called by CreateAndStartWorker

};
    
//                  G4WorkerRunManager* rm = GetRunManager();
#endif //G4UserWorkerInitialization_hh

