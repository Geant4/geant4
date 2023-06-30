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

  public:
    static G4WorkerRunManager* GetWorkerRunManager();
    static G4WorkerRunManagerKernel* GetWorkerRunManagerKernel();

    G4WorkerRunManager();
    ~G4WorkerRunManager() override;

    void InitializeGeometry() override;
    void RunInitialization() override;
    void DoEventLoop(G4int n_event, const char* macroFile = nullptr, G4int n_select = -1) override;
    void ProcessOneEvent(G4int i_event) override;
    G4Event* GenerateEvent(G4int i_event) override;

    void RunTermination() override;
    void TerminateEventLoop() override;

    // This function is called by the thread function: it should
    // loop until some work is requested.
    virtual void DoWork();

    // Sets the worker context.
    inline void SetWorkerThread(G4WorkerThread* wc) { workerContext = wc; }

    void SetUserInitialization(G4VUserPhysicsList* userInit) override;
    void SetUserInitialization(G4VUserDetectorConstruction* userInit) override;
    void SetUserInitialization(G4VUserActionInitialization* userInit) override;
    void SetUserInitialization(G4UserWorkerInitialization* userInit) override;
    void SetUserInitialization(G4UserWorkerThreadInitialization* userInit) override;
    void SetUserAction(G4UserRunAction* userAction) override;
    void SetUserAction(G4VUserPrimaryGeneratorAction* userAction) override;
    void SetUserAction(G4UserEventAction* userAction) override;
    void SetUserAction(G4UserStackingAction* userAction) override;
    void SetUserAction(G4UserTrackingAction* userAction) override;
    void SetUserAction(G4UserSteppingAction* userAction) override;

    void RestoreRndmEachEvent(G4bool flag) override { readStatusFromFile = flag; }

  protected:
    void ConstructScoringWorlds() override;
    void StoreRNGStatus(const G4String& filenamePrefix) override;
    void rndmSaveThisRun() override;
    void rndmSaveThisEvent() override;

    // This method will merge (reduce) the results
    // of this run into the global run
    virtual void MergePartialResults();

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
    void SetupDefaultRNGEngine();

  private:
    std::unique_ptr<ProfilerConfig> workerRunProfiler;
};

#endif  // G4WorkerRunManager_hh
