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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#pragma once

#include "G4ITReactionChange.hh"
#include "G4ITType.hh"
#include "G4ITModelHandler.hh"
#include "G4ITStepStatus.hh"
#include <vector>

class G4VITTimeStepComputer;
class G4VITReactionProcess;
class G4ITModelHandler;
class G4ITReactionSet;
class G4UserTimeStepAction;
class G4ITTrackingManager;
class G4ITTrackHolder;

/**
 * The G4ITModelProcessor will call the two processes defined in G4VITModel.
 * This processes act at the beginning and end of each step.
 * The first one, the TimeStepper will calculate a time step to propagate all
 * the track and eventually it can return some tracks that can likely react
 * at the end of the step.
 * The second one, the ReactionProcess will make the tracks reacting.
 * \deprecated This class will be removed
 */
class G4ITModelProcessor
{
public:
    G4ITModelProcessor();
    G4ITModelProcessor(const G4ITModelProcessor& other) = delete;
    G4ITModelProcessor& operator=(const G4ITModelProcessor& other) = delete;
    virtual ~G4ITModelProcessor();

    void SetModelHandler(G4ITModelHandler*);
    void SetTrackingManager(G4ITTrackingManager* trackingManager);

    void Initialize();

    void RegisterModel(double time, G4VITStepModel*);

    /** Restore the original state. This method should be called only by G4Scheduler */
    void CleanProcessor();

    G4double CalculateMinTimeStep(G4double currentGlobalTime,
                                  G4double definedMinTimeStep);

    void ComputeTrackReaction(G4ITStepStatus fITStepStatus,
                              G4double fGlobalTime,
                              G4double currentTimeStep,
                              G4double previousTimeStep,
                              G4bool reachedUserTimeLimit,
                              G4double fTimeTolerance,
                              G4UserTimeStepAction* fpUserTimeStepAction,
                              G4int fVerbose);

    void InitializeStepper(G4double currentGlobalTime,
                           G4double userMinTime);

    bool GetComputeTimeStep() const;

public:
    void CalculateTimeStep(const G4Track*, G4double userMinTimeStep);

    void DoCalculateStep();

    void FindReaction(G4ITReactionSet* pReactionSet,
                      double currentStepTime,
                      double previousStepTime,
                      bool reachedUserStepTimeLimit);

    const G4Track* GetTrack() const;

protected:
    void SetTrack(const G4Track*);
    void ExtractTimeStepperData();

    G4double fTSTimeStep;
    G4ITReactionSet* fReactionSet;
    G4ITTrackingManager* fpTrackingManager;
    G4ITTrackHolder* fpTrackContainer;

    G4bool fInitialized;
    G4ITModelHandler* fpModelHandler;

    const G4Track* fpTrack;
    G4double fUserMinTimeStep;

    std::vector<G4VITStepModel*> fActiveModels;
    G4VITStepModel* fpActiveModelWithMinTimeStep;

    std::vector<std::unique_ptr<G4ITReactionChange>> fReactionInfo;

    bool fComputeTimeStep;
    bool fComputeReaction;
};

