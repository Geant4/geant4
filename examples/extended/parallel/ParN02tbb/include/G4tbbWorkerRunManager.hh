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
//
// class description:
//
//      This is a class for run control in GEANT4 for TBB-type multi-
// threaded runs.
// It extends G4tbbWorkerRunManager 
// 
// Note: the User should never initialize instances of this class, that are 
//  usually handled by the TBB extension classes.
// There exists one instance of this class for each worker in a TBB application.
// 

#ifndef G4tbbWorkerRunManager_h
#define G4tbbWorkerRunManager_h 1

#include "G4WorkerRunManager.hh"
// #include "G4VScoringMesh.hh"

// #include <tbb/task.h>
#include <tbb/concurrent_queue.h>

// #include <tbb/enumerable_thread_specific.h>

#include "G4WorkerThread.hh"

class G4tbbWorkerRunManager : public G4WorkerRunManager 
{

public:
    G4tbbWorkerRunManager();
    ~G4tbbWorkerRunManager();
    // virtual void InitializeGeometry(); 

    //Modified for worker behavior
    virtual void ProcessOneEvent(G4int i_event);
    G4bool  DoOneEvent(G4int i_event, G4String msg=G4String(""));
      // These deal with the needs of tbb

    virtual void DoEventLoop(G4int n_event,const char* macroFile=0,G4int n_select=-1);

    typedef tbb::concurrent_queue<G4tbbWorkerRunManager*> 
                                  G4tbbWorkerRunManagerInstancesType;
    static G4tbbWorkerRunManagerInstancesType& GetInstancesList(); 

    // virtual void RunTermination();
    static void DestroyWorkersAndCleanup(); 
    static unsigned int NumberOfWorkers();
protected:
    // virtual void ConstructScoringWorlds();
    // virtual void StoreRNGStatus(const G4String& filenamePrefix );

    // Do not forget to set the worker context - in the Worker Thread
    // void SetWorkerThread( G4WorkerThread* wc ) { workerContext = wc; }

    //Global static instance of a queue containing all instances of this objects
    static G4tbbWorkerRunManagerInstancesType instancesList;

private:
    // Worker thread context
    G4WorkerThread* fWorkerContext; 
};
#endif //G4tbbWorkerRunManager_h
