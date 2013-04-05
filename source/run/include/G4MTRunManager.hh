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
// Users initializes an instance of this class instead of G4RunManager
// to start a multi-threaded simulation. This class will instantiate and
// constrol one or more G4WorkerRunManager classes

#ifndef G4MTRunManager_h
#define G4MTRunManager_h 1

#include "G4RunManager.hh"
#include "G4Threading.hh"
#include <list>
#include <map>

class G4ScoringManager;
class G4VUserWorkerInitialization;
class G4WorkerRunManager;

//TODO: Split random number storage from this class

class G4MTRunManager : public G4RunManager {
public:
    G4MTRunManager();
    virtual ~G4MTRunManager();
    //New method
    void SetNumberThreads( G4int n ) { nworkers = n; }
    //Returns i-th seed
    static long GetSeed( G4int i );
    //Returns number of available seeds
    static G4int NumberOfAvailableSeeds() { return seedsnum; }
public:
    //Inherited methods to re-implement for MT case
    virtual void InitializeEventLoop(G4int n_event, const char* macroFile=0, G4int n_select=-1);
    //The following do not do anything for this runmanager
    virtual void TerminateOneEvent();
    virtual void ProcessOneEvent(G4int i_event);
    virtual void TerminateEventLoop();
    virtual void ConstructScoringWorlds();
protected:
    //Initialize the seeds list, if derived class does not implement this method
    //A default generation will be used (nevents*2 random seeds)
    virtual void InitializeSeeds( G4int /*nevts*/) { };
    //Adds one seed to the list of seeds
    void AddOneSeed( long seed ) { seeds[seedsnum++] = seed; }
    void InitializeSeedsQueue( G4int ns ) { seeds = new long[ns]; }
    virtual void PrepareCommandsStack();
public:
    std::vector<G4String> GetCommandStack();
    //This method is invoked just before spawning the threads to
    //collect from UI managere the list of commands that threads
    //will execute.
private:
    G4int nworkers;
    //List of workers (i.e. thread)
    typedef std::list<G4Thread*> G4ThreadsList;
    G4ThreadsList threads;
    //Array of random seeds number
    static long* seeds;
    static G4int seedsnum;
    //List of workers run managers
    typedef std::list<G4WorkerRunManager*> G4WorkerRunManagerList;
    G4WorkerRunManagerList workersRM;
    //List of all workers run managers
    void DestroyWorkers();
    //Empty the workersList
    std::vector<G4String> uiCmdsForWorkers;
    //List of UI commands for workers.
    CLHEP::HepRandomEngine* masterRNGEngine;
    //Pointer to the mastet thread random engine
public:
    typedef std::map<G4int,G4VPhysicalVolume*> masterWorlds_t;
    static G4ScoringManager* GetMasterScoringManager() { return masterScM; }
    static masterWorlds_t& GetMasterWorlds() { return masterWorlds; }
    static void addWorld( G4int counter, G4VPhysicalVolume* w) { masterWorlds.insert( std::make_pair<G4int,G4VPhysicalVolume*>(counter,w) ); }
    void AddWorkerRunManager( G4WorkerRunManager* workerRM );
    //Should be called by Workers
    const CLHEP::HepRandomEngine* getMasterRandomEngine() const { return masterRNGEngine; }
private:
    //Handling of master thread scoring worlds, access to it is needed by workers
    static G4ScoringManager* masterScM;
    static masterWorlds_t masterWorlds;
    //Singleton implementing master thread behavior
    static G4MTRunManager* masterRM;
public: // with description
    static G4MTRunManager* GetMasterRunManager();
    // Returns the singleton instance of the run manager common to all threads implementing the master behavior
    static G4RunManagerKernel* GetMasterRunManagerKernel();
    // Returns the singleton instance of the run manager kernel common to all threads
    
    virtual void SetUserInitialization(G4VUserWorkerInitialization* userInit);
    virtual void SetUserInitialization(G4VUserPhysicsList* userPL);
    virtual void SetUserInitialization(G4VUserDetectorConstruction* userDC);
    virtual void SetUserAction(G4UserRunAction* userAction);
    virtual void SetUserAction(G4VUserPrimaryGeneratorAction* userAction);
    virtual void SetUserAction(G4UserEventAction* userAction);
    virtual void SetUserAction(G4UserStackingAction* userAction);
    virtual void SetUserAction(G4UserTrackingAction* userAction);
    virtual void SetUserAction(G4UserSteppingAction* userAction);

public: 
    // To be invoked solely from G4WorkerRunManager to merge the results
    void MergeScores(const G4ScoringManager* localScoringManager);
    void MergeRun(const G4Run* localRun);
};

#endif //G4MTRunManager_h
