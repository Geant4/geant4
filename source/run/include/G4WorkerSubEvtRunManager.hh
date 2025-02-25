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
// Author: Makoto Asai (JLAB) - 2024
//
// class description:
//
//  This is a class for run control in GEANT4 for multi-threaded runs
//  of worker thread running in sub-event parallel mode.
//  It extends G4RunManager re-implementing multi-threaded behavior in
//  key methods. See documentation for G4RunManager. User should never
//  initialize instances of this class, that are usually handled by
//  G4MTRunManager. There exists one instance of this class for each
//  worker in a MT application.

#ifndef G4WorkerSubEvtRunManager_h
#define G4WorkerSubEvtRunManager_h 1

#include "G4WorkerTaskRunManager.hh"

class G4WorkerThread;
class G4WorkerTaskRunManagerKernel;

// N.B. G4WorkerTaskRunManagerKernel is still used by this G4WorkerSubEvtRunManager
//  without any modification. -- future work for phase-II developments
using G4WorkerSubEvtRunManagerKernel = G4WorkerTaskRunManagerKernel;

class G4WorkerSubEvtRunManager : public G4WorkerTaskRunManager 
{
  public:
    static G4WorkerSubEvtRunManager* GetWorkerRunManager();
    static G4WorkerSubEvtRunManagerKernel* GetWorkerRunManagerKernel();
    G4WorkerSubEvtRunManager(G4int subEventType = 0);
    //G4WorkerSubEvtRunManager() = delete;

    // Modified for worker behavior
    void RunInitialization() override;
    void DoEventLoop(G4int n_event, const char* macroFile = nullptr, G4int n_select = -1) override;
    void ProcessOneEvent(G4int i_event) override;
    G4Event* GenerateEvent(G4int i_event) override;
    void RunTermination() override;
    void TerminateEventLoop() override;
    void DoWork() override;
    void RestoreRndmEachEvent(G4bool flag) override { readStatusFromFile = flag; }

    void DoCleanup() override;
    void ProcessUI() override;

    void SetUserInitialization(G4VUserPhysicsList* userPL) override;
    void SetUserInitialization(G4VUserDetectorConstruction* userDC) override;
    void SetUserInitialization(G4UserWorkerInitialization* userInit) override;
    void SetUserInitialization(G4UserWorkerThreadInitialization* userInit) override;
    void SetUserInitialization(G4VUserActionInitialization* userInit) override;
    void SetUserAction(G4UserRunAction* userAction) override;
    void SetUserAction(G4VUserPrimaryGeneratorAction* userAction) override;
    void SetUserAction(G4UserEventAction* userAction) override;
    void SetUserAction(G4UserStackingAction* userAction) override;
    void SetUserAction(G4UserTrackingAction* userAction) override;
    void SetUserAction(G4UserSteppingAction* userAction) override;

  protected:
    void StoreRNGStatus(const G4String& filenamePrefix) override;
    void SetupDefaultRNGEngine() override;

  private:
    G4int fSubEventType = -1;

  public:
    G4int GetSubEventType() const override
    { return fSubEventType; }
    void SetSubEventType(G4int) override;
};

#endif  // G4WorkerSubEvtRunManager_h
