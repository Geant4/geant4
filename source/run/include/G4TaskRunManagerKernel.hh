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
// Author: Jonathan Madsen (May 28st 2020)
//
// class description:
//
//     This is a class for mandatory control of GEANT4 kernel.
//     This class implements Worker behavior in a MT application.
//
//     This class is constructed by G4TaskRunManager. If a user uses his/her own
//     class instead of G4TaskRunManager, this class must be instantiated by
//     him/herself at the very beginning of the application and must be deleted
//     at the very end of the application. Also, following methods must be
//     invoked in the proper order.
//       DefineWorldVolume
//       InitializePhysics
//       RunInitialization
//       RunTermination
//
//     User must provide his/her own classes derived from the following
//     abstract class and register it to the RunManagerKernel.
//        G4VUserPhysicsList - Particle types, Processes and Cuts
//
//     G4TaskRunManagerKernel does not have any eveny loop. Handling of events
//     is managed by G4RunManager.
//
//     This class re-implements only the method that require special treatment
//     to implement worker behavior

#ifndef G4TaskRunManagerKernel_hh
#define G4TaskRunManagerKernel_hh 1

#include "rundefs.hh"
#include "G4RunManagerKernel.hh"
#include "G4TaskRunManager.hh"
#include "G4Threading.hh"

class G4WorkerThread;
class G4WorkerTaskRunManager;
#include <vector>

class G4TaskRunManagerKernel : public G4RunManagerKernel
{
 public:
  G4TaskRunManagerKernel();
  virtual ~G4TaskRunManagerKernel();

 public:  // with descroption
  static G4WorkerThread* GetWorkerThread();
  // G4ThreadPool tasks
  static void InitializeWorker();
  static void ExecuteWorkerInit();
  static void ExecuteWorkerTask();
  static void TerminateWorkerRunEventLoop();
  static void TerminateWorker();
  static void TerminateWorkerRunEventLoop(G4WorkerTaskRunManager*);
  static void TerminateWorker(G4WorkerTaskRunManager*);

  static std::vector<G4String>& InitCommandStack();
  // Fill decay tables with particle definition pointers of
  // decay products. This method has to be invoked by
  // MTRunManager before event loop starts on workers.
  void SetUpDecayChannels();
  // This method should be invoked by G4TaskRunManager
  void BroadcastAbortRun(G4bool softAbort);

 protected:
  void SetupShadowProcess() const;
  G4RUN_DLL static std::vector<G4String> initCmdStack;
};

#endif  // G4TaskRunManagerKernel_hh
