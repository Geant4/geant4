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

#ifndef G4RTWorkerInitialization_hh
#define G4RTWorkerInitialization_hh

class G4UserRunAction;
class G4VUserPrimaryGeneratorAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4RTRunAction;
class G4RTPrimaryGeneratorAction;
class G4RTTrackingAction;
class G4RTSteppingAction;

#include "G4Threading.hh"
#include "G4UserWorkerInitialization.hh"

class G4RTWorkerInitialization  : public G4UserWorkerInitialization
{
public: // with description
    G4RTWorkerInitialization();
    virtual ~G4RTWorkerInitialization();

    virtual void WorkerRunStart() const;
    // This method is called before an event loop. Geometry and physics have
    // already been set up for the thread. All threads are synchronized and
    // ready to start the local event loop. This situation is identical to
    // "Idle" state in the sequential mode.

    virtual void WorkerRunEnd() const;
    // This method is called for each thread, when the local event loop has
    // finished but before the synchronization over threads.

private:
    static G4ThreadLocal const G4UserRunAction * theUserRunAction;
    static G4ThreadLocal const G4VUserPrimaryGeneratorAction * theUserPrimaryGeneratorAction;
    static G4ThreadLocal const G4UserEventAction * theUserEventAction;
    static G4ThreadLocal const G4UserStackingAction * theUserStackingAction;
    static G4ThreadLocal const G4UserTrackingAction * theUserTrackingAction;
    static G4ThreadLocal const G4UserSteppingAction * theUserSteppingAction;

    static G4ThreadLocal G4RTRunAction * theRTRunAction;
    static G4ThreadLocal G4RTPrimaryGeneratorAction * theRTPrimaryGeneratorAction;
    static G4ThreadLocal G4RTTrackingAction * theRTTrackingAction;
    static G4ThreadLocal G4RTSteppingAction * theRTSteppingAction;

};
    
#endif //G4RTWorkerInitialization_hh

