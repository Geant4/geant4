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

#include "G4ITStepProcessorState_Lock.hh"
#include "G4ITLeadingTracks.hh"

class G4ITNavigator;
class G4ParticleDefinition;
class G4ITTrackingManager;
class G4IT;
class G4TrackingInformation;
class G4ITTransportation;
class G4VITProcess;
class G4VITSteppingVerbose;
class G4ITTrackHolder;
typedef class std::vector<int, std::allocator<int> > G4SelectedAtRestDoItVector;
typedef class std::vector<int, std::allocator<int> > G4SelectedAlongStepDoItVector;
typedef class std::vector<int, std::allocator<int> > G4SelectedPostStepDoItVector;

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
  // * Max number of processes
  std::size_t MAXofAtRestLoops;
  std::size_t MAXofAlongStepLoops;
  std::size_t MAXofPostStepLoops;
  // Maximum number of processes for each type of process
  // These depend on the G4ParticleDefinition, so on the track

  // * Transportation process
  G4ITTransportation* fpTransportation;
};

//________________________________________________
//
//          Members proper to a track
//________________________________________________
class G4ITStepProcessorState : public G4ITStepProcessorState_Lock
{
public:
  G4ITStepProcessorState();
  virtual ~G4ITStepProcessorState();

  G4ITStepProcessorState(const G4ITStepProcessorState&);
  G4ITStepProcessorState& operator=(const G4ITStepProcessorState&);

  // * Max Number of Process
  G4SelectedAtRestDoItVector fSelectedAtRestDoItVector;
  G4SelectedPostStepDoItVector fSelectedPostStepDoItVector;

  G4double fPhysicalStep;
  G4double fPreviousStepSize;
  G4double fSafety;

  G4StepStatus fStepStatus;

  // * Safety
  G4double fProposedSafety;
  // This keeps the minimum safety value proposed by AlongStepGPILs.
  G4ThreeVector fEndpointSafOrigin;
  G4double fEndpointSafety;
  // To get the true safety value at the PostStepPoint, you have
  // to subtract the distance to 'endpointSafOrigin' from this value.

  G4TouchableHandle fTouchableHandle;
};

/**
 * Its role is the same as G4StepManager :
 * - Find the minimum physical length and corresponding time step
 * - Step one track BUT on a given time step.
 */

class G4ITStepProcessor
{
friend class G4Scheduler;
public:
  G4ITStepProcessor();
  virtual ~G4ITStepProcessor();

  inline void SetPreviousStepTime(G4double);

  inline G4Track* GetTrack()
  {
    return fpTrack;
  }
  inline G4Step* GetStep()
  {
    return fpStep;
  }
  inline const G4Step* GetStep() const
  {
    return fpStep;
  }
  inline void SetStep(G4Step* val)
  {
    fpStep = val;
  }

  inline G4TrackVector* GetSecondaries() const
  {
    return fpSecondary;
  }
  inline void SetTrackingManager(G4ITTrackingManager* trackMan)
  {
    fpTrackingManager = trackMan;
  }
  inline G4ITTrackingManager* GetTrackingManager()
  {
    return fpTrackingManager;
  }

  //___________________________________

  virtual void Initialize();
  void ForceReInitialization();

  void ResetLeadingTracks();
  void PrepareLeadingTracks();

  //___________________________________
  G4double ComputeInteractionLength(double previousTimeStep);
  void DefinePhysicalStepLength(G4Track*);
  G4double GetILTimeStep()
  {
    return fILTimeStep;
  }

  //___________________________________
  // DoIt
  void DoIt(double timeStep);
  void ExtractDoItData();
  void Stepping(G4Track*, const double&);
  void FindTransportationStep();
  //___________________________________

  inline double GetInteractionTime();
  inline const G4Track* GetTrack() const;
  inline void CleanProcessor();

  std::size_t GetAtRestDoItProcTriggered() const
  {
    return fAtRestDoItProcTriggered;
  }

  G4GPILSelection GetGPILSelection() const
  {
    return fGPILSelection;
  }

  G4int GetN2ndariesAlongStepDoIt() const
  {
    return fN2ndariesAlongStepDoIt;
  }

  G4int GetN2ndariesAtRestDoIt() const
  {
    return fN2ndariesAtRestDoIt;
  }

  G4int GetN2ndariesPostStepDoIt() const
  {
    return fN2ndariesPostStepDoIt;
  }

  const G4VITProcess* GetCurrentProcess() const
  {
    return fpCurrentProcess;
  }

  G4double GetPhysIntLength() const
  {
    return fPhysIntLength;
  }

  std::size_t GetPostStepAtTimeDoItProcTriggered() const
  {
    return fPostStepAtTimeDoItProcTriggered;
  }

  std::size_t GetPostStepDoItProcTriggered() const
  {
    return fPostStepDoItProcTriggered;
  }

  const ProcessGeneralInfo* GetCurrentProcessInfo() const
  {
    return fpProcessInfo;
  }

  const G4ITStepProcessorState* GetProcessorState() const
  {
    return fpState;
  }

  const G4VParticleChange* GetParticleChange() const
  {
    return fpParticleChange;
  }

  const G4VPhysicalVolume* GetCurrentVolume() const
  {
    return fpCurrentVolume;
  }

  G4ForceCondition GetCondition() const
  {
    return fCondition;
  }

protected:

  void ExtractILData();

  void SetupGeneralProcessInfo(G4ParticleDefinition*, G4ProcessManager*);
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
  void PushSecondaries();

  // void CalculateStep();
  // void DoCalculateStep();

  // void CloneProcesses();
  void ActiveOnlyITProcess();
  void ActiveOnlyITProcess(G4ProcessManager*);

  void DealWithSecondaries(G4int&);
  void InvokeAtRestDoItProcs();
  void InvokeAlongStepDoItProcs();
  void InvokePostStepDoItProcs();
  void InvokePSDIP(std::size_t); //
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

  G4double kCarTolerance;
  // Cached geometrical tolerance on surface

  G4ITNavigator* fpNavigator;
  G4int fStoreTrajectory;
  G4VITSteppingVerbose* fpVerbose;

  G4ITTrackHolder* fpTrackContainer;
  G4ITLeadingTracks fLeadingTracks;

  //________________________________________________
  //
  // Members used as temporaries (= not proper to a track)
  //________________________________________________

  G4double fTimeStep; // not proper to a track
  G4double fILTimeStep; // proper to a track ensemble

  G4double fPreviousTimeStep;
  G4TrackVector* fpSecondary; // get from fpStep at every configuration setup
  G4VParticleChange* fpParticleChange;

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
  std::size_t fAtRestDoItProcTriggered;
  std::size_t fPostStepDoItProcTriggered;
  std::size_t fPostStepAtTimeDoItProcTriggered;
  // Record the selected process

  G4ForceCondition fCondition;
  G4GPILSelection fGPILSelection;
  // Above three variables are for the method
  // DefinePhysicalStepLength(). To pass these information to
  // the method Verbose, they are kept at here. Need a more
  // elegant mechanism.

  G4double fPhysIntLength;
  // The minimum physical interaction length over all possible processes

  // * Sensitive detector
  //    G4SteppingControl StepControlFlag;
  //    G4VSensitiveDetector*   fpSensitive;

  G4VPhysicalVolume* fpCurrentVolume;
  // Get from fpStep or touchable, keep as member for user interface

  //________________________________________________
  //
  // Members related to ParticleDefinition and not
  // proper to a track
  //________________________________________________

  std::map<const G4ParticleDefinition*, ProcessGeneralInfo*> fProcessGeneralInfoMap;
  ProcessGeneralInfo* fpProcessInfo;
  G4ITTransportation* fpTransportation;

  //________________________________________________
  //
  // Members used for setting up the processor
  //________________________________________________

  G4Track* fpTrack; // Set track
  G4IT* fpITrack; // Set track
  G4TrackingInformation* fpTrackingInfo; // Set track

  G4ITStepProcessorState* fpState; // SetupMembers or InitDefineStep
  G4Step* fpStep; // Set track or InitDefineStep

  G4StepPoint* fpPreStepPoint; // SetupMembers
  G4StepPoint* fpPostStepPoint; // SetupMembers
};

//______________________________________________________________________________

inline void G4ITStepProcessor::SetPreviousStepTime(G4double previousTimeStep)
{
  fPreviousTimeStep = previousTimeStep;
}

//______________________________________________________________________________

inline const G4Track* G4ITStepProcessor::GetTrack() const
{
  return fpTrack;
}

//______________________________________________________________________________

inline G4double G4ITStepProcessor::CalculateSafety()
{
  return std::max(fpState->fEndpointSafety - (fpState->fEndpointSafOrigin
                      - fpPostStepPoint->GetPosition()).mag(),
                  kCarTolerance);
}

//______________________________________________________________________________

inline void G4ITStepProcessor::SetNavigator(G4ITNavigator *value)
{
  fpNavigator = value;
}

//______________________________________________________________________________

inline void G4ITStepProcessor::CleanProcessor()
{
  fTimeStep = DBL_MAX;
  fPhysIntLength = DBL_MAX;

  fpState = 0;
  fpTrack = 0;
  fpTrackingInfo = 0;
  fpITrack = 0;
  fpStep = 0;
  fpPreStepPoint = 0;
  fpPostStepPoint = 0;

  fpParticleChange = 0;

  fpCurrentVolume = 0;
  //    fpSensitive = 0;

  fpSecondary = 0;

  fpTransportation = 0;

  fpCurrentProcess= 0;
  fpProcessInfo = 0;

  fAtRestDoItProcTriggered = INT_MAX;
  fPostStepDoItProcTriggered = INT_MAX;
  fPostStepAtTimeDoItProcTriggered = INT_MAX;
  fGPILSelection = NotCandidateForSelection;
  fCondition = NotForced;
}

//______________________________________________________________________________

inline double G4ITStepProcessor::GetInteractionTime()
{
  return fTimeStep;
}

#endif // G4ITSTEPPROCESSOR_H
