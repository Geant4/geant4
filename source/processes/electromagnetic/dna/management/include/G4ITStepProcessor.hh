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
// $Id: G4ITStepProcessor.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4ITSTEPPROCESSOR_H
#define G4ITSTEPPROCESSOR_H

#include "G4ios.hh"                   // Include from 'system'
#include "globals.hh"                 // Include from 'global'
#include "Randomize.hh"               // Include from 'global'

#include "G4LogicalVolume.hh"         // Include from 'geometry'
#include "G4VPhysicalVolume.hh"       // Include from 'geometry'
#include "G4ProcessManager.hh"        // Include from 'piim'

#include "G4Track.hh"                 // Include from 'track'
#include "G4TrackVector.hh"           // Include from 'track'
#include "G4TrackStatus.hh"           // Include from 'track'
#include "G4StepStatus.hh"            // Include from 'track'
//#include "G4UserSteppingAction.hh"    // Include from 'tracking'
//#include "G4UserTrackingAction.hh"    // Include from 'tracking'
#include "G4Step.hh"                  // Include from 'track'
#include "G4StepPoint.hh"             // Include from 'track'
#include "G4TouchableHandle.hh"             // Include from 'geometry'
#include "G4TouchableHistoryHandle.hh"      // Include from 'geometry'


#include "G4TrackingInformation.hh"

//class G4Navigator;
class G4ITNavigator;
class G4ParticleDefinition ;
class G4ITTrackingManager;
class G4IT;
class G4TrackingInformation;
class G4ITTransportation;
class G4VITProcess;
typedef class std::vector<int, std::allocator<int> > G4SelectedAtRestDoItVector;
typedef class std::vector<int, std::allocator<int> > G4SelectedAlongStepDoItVector;
typedef class std::vector<int, std::allocator<int> > G4SelectedPostStepDoItVector;


/**
  * Its role is the same as G4StepManager :
  * - Find the minimum physical length and corresponding time step
  * - Step one track BUT on a given time step.
  */

class G4ITStepProcessor
{

public:
    G4ITStepProcessor();
    virtual ~G4ITStepProcessor();

    inline void SetPreviousStepTime(G4double);

    inline G4Track* GetTrack()                                            {return fpTrack;}
    inline G4Step* GetStep()                                              {return fpStep;}
    inline const G4Step* GetStep() const                                  {return fpStep;}
    inline void SetStep(G4Step* val)                                      {fpStep = val;}

    inline G4TrackVector* GetSecondaries()                         {return fpSecondary;}
    inline void SetTrackingManager(G4ITTrackingManager* trackMan)  {fpTrackingManager = trackMan;}
    inline G4ITTrackingManager* GetTrackingManager()               {return fpTrackingManager;}

    virtual void Initialize();
    void ForceReInitialization();

    void DefinePhysicalStepLength(G4Track*);
    void Stepping(G4Track*, const double&);
    void CalculateStep(G4Track*, const double&);
    void CalculateStep(G4Track*);

    void DoIt(G4Track*,double);

    void FindTransportationStep();
    void UpdateTrack(G4Track*);

    inline double GetInteractionTime();
    inline const G4Track* GetTrack() const ;
    inline void CleanProcessor();

protected:
    void SetupGeneralProcessInfo(G4ParticleDefinition*,G4ProcessManager*);
    void ClearProcessInfo();
    void SetTrack(G4Track*);

    void GetProcessInfo();

    void SetupMembers();
    void ResetSecondaries();
    void InitDefineStep();

    void SetInitialStep();

    void GetAtRestIL();
    void DoDefinePhysicalStepLength();
    void DoStepping();

    void CalculateStep();
    void DoCalculateStep();

    void CloneProcesses();
    void ActiveOnlyITProcess();
    void ActiveOnlyITProcess(G4ProcessManager*);

    void DealWithSecondaries(G4int&);
    void InvokeAtRestDoItProcs();
    void InvokeAlongStepDoItProcs();
    void InvokePostStepDoItProcs();
    void InvokePSDIP(size_t); //
    void InvokeTransportationProc();
    void SetNavigator(G4ITNavigator *value);
    G4double CalculateSafety();

    // Return the estimated safety value at the PostStepPoint
    void ApplyProductionCut(G4Track*);


    G4ITStepProcessor(const G4ITStepProcessor& other);
    G4ITStepProcessor& operator=(const G4ITStepProcessor& other);

private:
    //________________________________________________
    //
    //              General members
    //________________________________________________

    G4bool fInitialized;

    G4ITTrackingManager* fpTrackingManager;
    //  G4UserSteppingAction*   fpUserSteppingAction;

    G4double kCarTolerance;
    // Cached geometrical tolerance on surface

    G4ITNavigator*            fpNavigator;
//    G4Navigator*            fpNavigator;
    G4int                   fStoreTrajectory;
    G4int                   verboseLevel;

    //________________________________________________
    //
    // Members used as temporaries (= not proper to a track)
    //________________________________________________

    G4double                fTimeStep ; // not proper to a track
    G4double                fPreviousTimeStep;
    G4TrackVector*          fpSecondary ; // get from fpStep at every configuration setup
    G4VParticleChange*      fpParticleChange;

    G4VITProcess* fpCurrentProcess;
    // The pointer to the process of which DoIt or
    // GetPhysicalInteractionLength has been just executed

    // * Secondaries
    G4int fN2ndariesAtRestDoIt;
    G4int fN2ndariesAlongStepDoIt;
    G4int fN2ndariesPostStepDoIt;
    // These are the numbers of secondaries generated by the process
    // just executed.

    // * Process selection
    size_t fAtRestDoItProcTriggered;
    size_t fPostStepDoItProcTriggered;
    size_t fPostStepAtTimeDoItProcTriggered;
    // Record the selected process

    G4ForceCondition fCondition;
    G4GPILSelection  fGPILSelection;
    // Above three variables are for the method
    // DefinePhysicalStepLength(). To pass these information to
    // the method Verbose, they are kept at here. Need a more
    // elegant mechanism.

    G4double fPhysIntLength;
    // The minimum physical interaction length over all possible processes

    // * Sensitive detector
//    G4SteppingControl StepControlFlag;
//    G4VSensitiveDetector*   fpSensitive;

    G4VPhysicalVolume*      fpCurrentVolume; // Get from fpStep or touchable, keep as member for user interface

    //________________________________________________
    //
    // Members related to ParticleDefinition and not
    // proper to a track
    //________________________________________________
    struct ProcessGeneralInfo
    {
        G4ProcessVector* fpAtRestDoItVector;
        G4ProcessVector* fpAlongStepDoItVector;
        G4ProcessVector* fpPostStepDoItVector;

        G4ProcessVector* fpAtRestGetPhysIntVector;
        G4ProcessVector* fpAlongStepGetPhysIntVector;
        G4ProcessVector* fpPostStepGetPhysIntVector;
        //
        // Note: DoItVector has inverse order against GetPhysIntVector
        //       and SelectedPostStepDoItVector.
        //
        // * Max Number of Process
        size_t MAXofAtRestLoops;
        size_t MAXofAlongStepLoops;
        size_t MAXofPostStepLoops;
        // Maximum number of processes for each type of process
        // These depend on the G4ParticleDefinition, so on the track

        // * Transportation process
        G4ITTransportation* fpTransportation ;
    };

    std::map<const G4ParticleDefinition*, ProcessGeneralInfo*> fProcessGeneralInfoMap;
    ProcessGeneralInfo* fpProcessInfo;

    G4ITTransportation* fpTransportation ;

    //________________________________________________
    //
    //          Members proper to a track
    //________________________________________________
    class G4ITStepProcessorState : public G4ITStepProcessorState_Lock
    {
    public:
        G4ITStepProcessorState();
        virtual ~G4ITStepProcessorState();

        // * Max Number of Process
        G4SelectedAtRestDoItVector fSelectedAtRestDoItVector;
        G4SelectedPostStepDoItVector fSelectedPostStepDoItVector;

        G4double    fPhysicalStep;
        G4double    fPreviousStepSize;
        G4double    fSafety;

        G4StepStatus fStepStatus;

        // * Safety
        G4double proposedSafety;
        // This keeps the minimum safety value proposed by AlongStepGPILs.
        G4ThreeVector endpointSafOrigin;
        G4double endpointSafety;
        // To get the true safety value at the PostStepPoint, you have
        // to subtract the distance to 'endpointSafOrigin' from this value.

        G4TouchableHandle fTouchableHandle;
    private :
        G4ITStepProcessorState(const G4ITStepProcessorState&);
        G4ITStepProcessorState&  operator=(const G4ITStepProcessorState&);
    };

    //________________________________________________
    //
    // Members used for configurating the processor
    //________________________________________________

    G4Track*                fpTrack; // Set track
    G4IT*                   fpITrack ; // Set track
    G4TrackingInformation*  fpTrackingInfo ; // Set track

    G4ITStepProcessorState* fpState; // SetupMembers or InitDefineStep
    G4Step*                 fpStep; // Set track or InitDefineStep

    G4StepPoint*            fpPreStepPoint; // SetupMembers
    G4StepPoint*            fpPostStepPoint; // SetupMembers
};

inline void G4ITStepProcessor::SetPreviousStepTime(G4double previousTimeStep)
{
    fPreviousTimeStep = previousTimeStep;
}

inline const G4Track* G4ITStepProcessor::GetTrack() const
{
    return fpTrack;
}

inline G4double G4ITStepProcessor::CalculateSafety()
{
    return std::max( fpState->endpointSafety -
                     (fpState->endpointSafOrigin - fpPostStepPoint->GetPosition()).mag(),
                     kCarTolerance );
}

inline void G4ITStepProcessor::SetNavigator(G4ITNavigator *value)
{
    fpNavigator = value;
}

inline void G4ITStepProcessor::CleanProcessor()
{
    fTimeStep = DBL_MAX ;
    fPhysIntLength = DBL_MAX;

    fpState = 0;
    fpTrack = 0;
    fpTrackingInfo = 0 ;
    fpITrack = 0;
    fpStep = 0;
    fpPreStepPoint = 0;
    fpPostStepPoint = 0;

    fpParticleChange = 0;

    fpCurrentVolume = 0;
//    fpSensitive = 0;

    fpSecondary = 0 ;

    fpTransportation = 0;

    fpCurrentProcess= 0;
    fpProcessInfo = 0;

    fAtRestDoItProcTriggered = INT_MAX;
    fPostStepDoItProcTriggered = INT_MAX;
    fPostStepAtTimeDoItProcTriggered = INT_MAX;
    fGPILSelection = NotCandidateForSelection ;
    fCondition = NotForced;
}

//______________________________________________________________________________
inline double G4ITStepProcessor::GetInteractionTime()
{
    return fTimeStep ;
}


#endif // G4ITSTEPPROCESSOR_H
