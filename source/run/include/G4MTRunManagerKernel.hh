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
////
//

// class description:
//
//     This is a class for mandatory control of GEANT4 kernel.
//     This class implements Worker behavior in a MT application.
//
//     This class is constructed by G4MTRunManager. If a user uses his/her own
//     class instead of G4MTRunManager, this class must be instantiated by
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
//     G4MTRunManagerKernel does not have any eveny loop. Handling of events
//     is managed by G4RunManager.
//
//     This class re-implements only the method that require special treatment
//     to implement worker behavior

#ifndef G4MTRunManagerKernel_h
#define G4MTRunManagerKernel_h 1

#include "G4RunManagerKernel.hh"
#include "G4Threading.hh"
#include "G4MTRunManager.hh"

class G4WorkerThread;
class G4WorkerRunManager;
#include <vector>

class G4MTRunManagerKernel : public G4RunManagerKernel {
public:
    G4MTRunManagerKernel();
    virtual ~G4MTRunManagerKernel();
protected:
    void SetupShadowProcess() const;

public: // with descroption
    // This static method is used to start a worker thread.
    // Virtual methods to be invoked from this methos are
    // defined in G4UserWorkerInitialization class.
    static void* StartThread(void* context);

//private:
//    static void ReinitializeGeometry();
private:
    static G4ThreadLocal G4WorkerThread* wThreadContext;
public:
    static G4WorkerThread* GetWorkerThread();

public: // with descroption
    // Fill decay tables with particle definition pointers of
    // decay products. This method has to be invoked by 
    // MTRunManager before event loop starts on workers.
    void SetUpDecayChannels();

private:
    static std::vector<G4WorkerRunManager*>* workerRMvector;

public:
    // This method should be invoked by G4MTRunManager
    void BroadcastAbortRun(G4bool softAbort);
};

#endif //G4MTRunManagerKernel_h
