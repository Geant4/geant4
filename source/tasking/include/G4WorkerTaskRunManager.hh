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
// Author: Jonathan Madsen (May 28st 2020)
//
// class description:
//
//  This is a class for run control in GEANT4 for multi-threaded runs
//  It extends G4RunManager re-implementing multi-threaded behavior in
//  key methods. See documentation for G4RunManager. User should never
//  initialize instances of this class, that are usually handled by
//  G4MTRunManager. There exists one instance of this class for each
//  worker in a MT application.

#ifndef G4WorkerTaskRunManager_h
#define G4WorkerTaskRunManager_h 1

#include "G4RNGHelper.hh"
#include "G4RunManager.hh"
#include "G4WorkerRunManager.hh"

class G4WorkerThread;
class G4WorkerTaskRunManagerKernel;

class G4WorkerTaskRunManager : public G4WorkerRunManager
{
 public:
  using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Run>;

 public:
  typedef std::vector<G4String> G4StrVector;

 public:
  static G4WorkerTaskRunManager* GetWorkerRunManager();
  static G4WorkerTaskRunManagerKernel* GetWorkerRunManagerKernel();
  G4WorkerTaskRunManager();

  // Modified for worker behavior
  virtual void RunInitialization() override;
  virtual void DoEventLoop(G4int n_event, const char* macroFile = nullptr,
                           G4int n_select = -1) override;
  virtual void ProcessOneEvent(G4int i_event) override;
  virtual G4Event* GenerateEvent(G4int i_event) override;
  virtual void RunTermination() override;
  virtual void TerminateEventLoop() override;
  virtual void DoWork() override;
  virtual void RestoreRndmEachEvent(G4bool flag) override
  {
    readStatusFromFile = flag;
  }

  virtual void DoCleanup();
  virtual void ProcessUI();
  G4WorkerThread* GetWorkerThread() const { return workerContext; }
  G4StrVector GetCommandStack() const { return processedCommandStack; }

 protected:
  virtual void StoreRNGStatus(const G4String& filenamePrefix) override;

 private:
  void SetupDefaultRNGEngine();

 private:
  G4StrVector processedCommandStack;

 private:
  std::unique_ptr<ProfilerConfig> workerRunProfiler;
};

#endif  // G4WorkerTaskRunManager_h
