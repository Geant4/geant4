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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4ITStepManager_h
#define G4ITStepManager_h

#include "G4TrackList.hh"
#include <globals.hh>
#include <vector>
#include <map>
#include "G4ITModelHandler.hh"

class G4ITTrackingManager;
class G4ITModelProcessor;
class G4ITStepProcessor;
class G4Track ;
class G4UserReactionAction;

class G4ITStepManager
{

public :
    static G4ITStepManager* Instance();

    ~G4ITStepManager();

    void Initialize() ;
    void Reset();
    void Process() ;
    void PushTrack(G4Track*);
    void ClearList();
    void SetEndTime(const double);

    inline G4ITModelHandler* GetModelHandler();
    inline void SetTimeSteps(std::map<double,double>*);
    inline G4int GetNbSteps();
    inline G4double GetEndTime();
    inline G4double GetMinTimeStep();
    inline G4double GetGlobalTime();
    inline void SetUserITAction(G4UserReactionAction* /*userITAction*/);
    inline G4UserReactionAction* GetUserITAction();

    inline void SetVerbose(int);
    inline int GetVerbose();

protected:
    void DoProcess();
    void SynchronizeTracks();
    void Stepping();

    void FindUserPreDefinedTimeStep();
    void CalculateMinStep() ;
    void ComputeInteractionLength();
    void DoIt();
    void ComputeInteractionBetweenTracks();
    void MergeSecondariesWithMainList();

    void PushSecondaries(G4ITStepProcessor*);
    void _PushTrack(G4Track*);
    void PushDelayed(G4Track*, const G4double&);

    void EndTracking(G4Track*);
    void KillTracks();
    void EndTracking();

    void ExtractTimeStepperData(G4ITModelProcessor*);
    void ExtractILData(G4ITStepProcessor*);
    void ExtractDoItData(G4ITStepProcessor*);

private:
    G4ITStepManager();

    bool fInitialized;
    bool fRunning ;

    // Time members
    double fGlobalTime;
    double fTmpGlobalTime ;
    double fEndTime ;       // fEndTime : stores the biggest global time steps here (for synchronizing tracks)
    double fTmpEndTime ;    // fTmpEndTime : Stores the biggest end time here (for synchronizing tracks)

    // Flags
    bool fComputeTimeStep;
    bool fComputeReaction;

    double fMinTimeStep ; // The selected minimum time step

    // User steps
    bool fUsePreDefinedTimeSteps ;
    std::map<double, double>* fpUserTimeSteps ; // One can give time steps in respect to the global time
    double fDefinedMinTimeStep  ; // selected user time step in respect to the global time
    bool fReachedUserTimeLimit ; // if fMinTimeStep == the user time step

    int fNbSteps ;

    std::map<G4Track*, G4TrackVectorHandle> fReactingTracks ;
    std::vector<G4Track*> fLeadingTracks ;

    bool fInteractionStep ;
    // Flag : if the step is lead by the interaction with the matter and NOT by tracks' reaction

    G4TrackList*                    fpMainList ;
    G4TrackList                     fSecondaries ; //merge
    G4TrackList*                    fpWaitingList ; // Waiting queue of currentList
    std::map<double,G4TrackList* >  fDelayedList ;
    G4TrackList                     fToBeKilledList ;

    G4ITStepProcessor* fpStepProcessor ;
    G4ITStepProcessor* fpMasterStepProcessor ;

    G4ITModelProcessor* fpModelProcessor;
    G4ITModelProcessor* fpMasterModelProcessor;

    G4ITModelHandler* fpModelHandler;

    G4ITTrackingManager* fpTrackingManager;
    G4UserReactionAction* fpUserITAction;

    int fNbTracks ;

    int fVerbose;

    static G4ITStepManager* fgStepManager ;
};

inline G4ITModelHandler* G4ITStepManager::GetModelHandler() {return fpModelHandler;}

inline void G4ITStepManager::SetEndTime(const double __endtime) { fEndTime = __endtime ;}

inline void G4ITStepManager::SetTimeSteps(std::map<double,double>* steps)
{
    fUsePreDefinedTimeSteps = true;
    fpUserTimeSteps = steps ;
}

inline G4int G4ITStepManager::GetNbSteps()
{
    return fNbSteps;
}

inline G4double G4ITStepManager::GetEndTime()
{
    return fEndTime;
} 

inline G4double G4ITStepManager::GetMinTimeStep()
{
    return fMinTimeStep ;
}


inline G4double G4ITStepManager::GetGlobalTime()
{
    return fGlobalTime ;
}

inline void G4ITStepManager::SetUserITAction(G4UserReactionAction* userITAction)
{
    fpUserITAction = userITAction;
}

inline G4UserReactionAction* G4ITStepManager::GetUserITAction()
{
    return fpUserITAction;
}

inline void G4ITStepManager::SetVerbose(int verbose)
{
    fVerbose = verbose;
}

inline int G4ITStepManager::GetVerbose()
{
    return fVerbose ;
}

#endif
