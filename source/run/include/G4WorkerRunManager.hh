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
// G4WorkerRunManager
//
// Class description:
//
// This is a class for run control in GEANT4 for multi-threading.
// It extends G4RunManager re-implementing multi-threaded behavior in
// key methods. See documentation for G4RunManager.
// User should never initialize instances of this class, that are usually
// handled by G4MTRunManager. There exists one instance of this class for
// each worker in a MT application.

// Original authors: X.Dong, A.Dotti - 2013
// --------------------------------------------------------------------
#ifndef G4WorkerRunManager_hh
#define G4WorkerRunManager_hh 1

#include "G4RNGHelper.hh"
#include "G4RunManager.hh"

class G4WorkerThread;
class G4WorkerRunManagerKernel;

class G4WorkerRunManager : public G4RunManager
{
  public:

    using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Run>;

    static G4WorkerRunManager* GetWorkerRunManager();
    static G4WorkerRunManagerKernel* GetWorkerRunManagerKernel();

    G4WorkerRunManager();
   ~G4WorkerRunManager();

    virtual void InitializeGeometry();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event, const char* macroFile = 0,
                             G4int n_select = -1);
    virtual void ProcessOneEvent(G4int i_event);
    virtual G4Event* GenerateEvent(G4int i_event);

    virtual void RunTermination();
    virtual void TerminateEventLoop();

    virtual void DoWork();
      // This function is called by the thread function: it should
      // loop until some work is requested.

    inline void SetWorkerThread(G4WorkerThread* wc) { workerContext = wc; }
      // Sets the worker context.

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

    virtual void RestoreRndmEachEvent(G4bool flag) { readStatusFromFile=flag; }

  protected:

    virtual void ConstructScoringWorlds();
    virtual void StoreRNGStatus(const G4String& filenamePrefix);
    virtual void rndmSaveThisRun();
    virtual void rndmSaveThisEvent();
    virtual void MergePartialResults();
      // This method will merge (reduce) the results
      // of this run into the global run

  private:

    void SetupDefaultRNGEngine();

  protected:

    G4WorkerThread* workerContext = nullptr;
  #ifdef G4MULTITHREADED
    G4bool visIsSetUp = false;
  #endif

    G4bool eventLoopOnGoing = false;
    G4bool runIsSeeded = false;
    G4int nevModulo = -1;
    G4int currEvID = -1;
    G4int luxury = -1;
    G4SeedsQueue seedsQueue;
    G4bool readStatusFromFile = false;

 private:

  std::unique_ptr<ProfilerConfig> workerRunProfiler;
};

#endif  // G4WorkerRunManager_hh
