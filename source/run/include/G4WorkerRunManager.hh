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
//

// class description:
//
//      This is a class for run control in GEANT4 for multi-threaded runs
// It extends G4RunManager re-implementing multi-threaded behavior in
// key methods. See documentation for G4RunManager
// User should never initialize instances of this class, that are usually handled
// by G4MTRunManager.
// There exists one instance of this class for each worker in a MT application.

#ifndef G4WorkerRunManager_h
#define G4WorkerRunManager_h 1
#include "G4RNGHelper.hh"
#include "G4RunManager.hh"

class G4WorkerThread;
class G4WorkerRunManagerKernel;

class G4WorkerRunManager : public G4RunManager {
public:
    static G4WorkerRunManager* GetWorkerRunManager();
    static G4WorkerRunManagerKernel* GetWorkerRunManagerKernel();
    G4WorkerRunManager();
    ~G4WorkerRunManager();
    //Modified for worker behavior
    ////////virtual void BeamOn(G4int n_event,const char* macroFile=0,G4int n_select=-1);
    virtual void InitializeGeometry();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event,const char* macroFile=0,G4int n_select=-1);
    virtual void ProcessOneEvent(G4int i_event);
    virtual G4Event* GenerateEvent(G4int i_event);
    //G4int NewCommands( const std::vector<G4String>& newCmdsToExecute , G4String& currentCmd );
    //Called by the MTRunManager when new UI commands are to be executed.
    //It is not assumed method is not thread-safe: i.e. should be called sequentially
    //Returns 0 if commands are executed corrected, otherwise returns error code (see G4UImanager::ApplyCommand)
    //In case of error currentCmd is set to the command that gave the problem
    virtual void RunTermination();
    virtual void TerminateEventLoop();

    //This function is called by the thread function: it should loop until some
    //work is requested
    virtual void DoWork();
protected:
    virtual void ConstructScoringWorlds();
    virtual void StoreRNGStatus(const G4String& filenamePrefix );
    virtual void MergePartialResults();
    //This method will merge (reduce) the results of this run into the
    //global run
public:
    //! Sets the worker context
        void SetWorkerThread( G4WorkerThread* wc ) { workerContext = wc; }
private:
    G4WorkerThread* workerContext;
    void SetupDefaultRNGEngine();

public:
    virtual void SetUserInitialization(G4VUserPhysicsList* userInit);
    virtual void SetUserInitialization(G4VUserDetectorConstruction* userInit);
    virtual void SetUserInitialization(G4VUserActionInitialization* userInit);
    virtual void SetUserInitialization(G4UserWorkerInitialization* userInit);
    virtual void SetUserInitialization(G4UserWorkerThreadInitialization* userInit);
    virtual void SetUserAction(G4UserRunAction* userAction);
    virtual void SetUserAction(G4VUserPrimaryGeneratorAction* userAction);
    virtual void SetUserAction(G4UserEventAction* userAction);
    virtual void SetUserAction(G4UserStackingAction* userAction);
    virtual void SetUserAction(G4UserTrackingAction* userAction);
    virtual void SetUserAction(G4UserSteppingAction* userAction);

protected:
    G4bool eventLoopOnGoing;
    G4bool runIsSeeded;
    G4int nevModulo;
    G4int currEvID;
    G4SeedsQueue seedsQueue;
    G4bool readStatusFromFile;

public:
    virtual void RestoreRndmEachEvent(G4bool flag) { readStatusFromFile = flag; }
#ifdef G4MULTITHREADED
private:
    G4bool visIsSetUp;
#endif
};

#endif //G4WorkerRunManager_h
