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
// $Id: G4ITStepManager.hh 60427 2012-07-11 16:34:35Z matkara $
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

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

#ifndef G4Scheduler_h
#define G4Scheduler_h

#include <G4VScheduler.hh>
#include <vector>
#include <map>
#include <memory>

#include "globals.hh"

#include "G4ITModelHandler.hh"
#include "G4ITStepStatus.hh"
#include "G4ITTrackHolder.hh"
#include "G4VStateDependent.hh"
#include "G4ITReaction.hh"

class G4ITTrackingManager;
class G4ITModelProcessor;
class G4ITStepProcessor;
class G4Track;
class G4UserTimeStepAction;
class G4SchedulerMessenger;
class G4ITTrackingInteractivity;
class G4ITGun;

#ifndef compTrackPerID__
#define compTrackPerID__
  struct compTrackPerID
  {
    bool operator()(G4Track* rhs, G4Track* lhs) const
    {
      return rhs->GetTrackID() < lhs->GetTrackID();
    }
  };
#endif


/**
 * G4ITStepManager enables to synchronize in time
 * the step of tracks.
 */
class G4Scheduler :
    public G4VScheduler,
    public G4VStateDependent
{
protected:
  virtual ~G4Scheduler();

public:
  static G4Scheduler* Instance();
  /** DeleteInstance should be used instead
   * of the destructor
   */
  static void DeleteInstance();
  virtual G4bool Notify(G4ApplicationState requestedState);

  virtual void RegisterModel(G4VITStepModel*, double);

  void Initialize();
  void ForceReinitialization();
  inline bool IsInitialized();
  inline bool IsRunning(){return fRunning;}
  void Reset();
  void Process();
  virtual void PushTrack(G4Track*);
  void ClearList();

  inline void SetGun(G4ITGun*);

  inline void Stop();
  void Clear();

  // To be called only in UserReactionAction::EndProcessing()
  // after fRunning flag has been turned off.
  // This is not done automatically before UserReactionAction::EndProcessing()
  // is called in case one would like to access some track information
  void EndTracking();

  void SetEndTime(const double);

  inline void SetTimeTolerance(double);
  // Two tracks below the time tolerance are supposed to be
  // in the same time slice
  inline double GetTimeTolerance() const;

  inline void SetMaxZeroTimeAllowed(int);
  inline int GetMaxZeroTimeAllowed() const;

  inline G4ITModelHandler* GetModelHandler();
  inline void SetTimeSteps(std::map<double, double>*);
  inline void AddTimeStep(double, double);
  inline void SetDefaultTimeStep(double);
  double GetLimitingTimeStep() const;
  inline G4int GetNbSteps() const;
  inline void SetMaxNbSteps(G4int);
  inline G4int GetMaxNbSteps() const;
  inline G4double GetStartTime() const;
  inline G4double GetEndTime() const;
  virtual inline G4double GetTimeStep() const;
  inline G4double GetPreviousTimeStep() const;
  inline G4double GetGlobalTime() const;
  inline void SetUserAction(G4UserTimeStepAction*);
  inline G4UserTimeStepAction* GetUserTimeStepAction() const;

  // To use with transportation only, no reactions
  inline void UseDefaultTimeSteps(G4bool);
  inline G4bool AreDefaultTimeStepsUsed();

  inline G4bool GetComputeTimeStepFlag() const;

  inline G4ITStepStatus GetStatus() const;

  inline void SetVerbose(int);
  // 1 : Reaction information
  // 2 : (1) + time step information
  // 3 : (2) + step info for individual tracks
  // 4 : (2) + trackList processing info + pushed and killed track info
  inline int GetVerbose() const;

  inline void WhyDoYouStop();

  void SetInteractivity(G4ITTrackingInteractivity*);
  inline G4ITTrackingInteractivity* GetInteractivity();

  virtual size_t GetNTracks();

  void GetCollisionType(G4String& interactionType);

protected:

  void DoProcess();
  void SynchronizeTracks();
  void Stepping();

  void FindUserPreDefinedTimeStep();
  void CalculateMinTimeStep();
  void ComputeInteractionLength();
  void DoIt();
  void ComputeTrackReaction();
//  void MergeSecondariesWithMainList();

  void PushSecondaries(G4ITStepProcessor*);
//  void _PushTrack(G4Track*);
  void PushDelayed(G4Track*, const G4double&);

  void EndTracking(G4Track*);
  void KillTracks();

  void ExtractTimeStepperData(G4ITModelProcessor*);
  void ExtractILData(G4ITStepProcessor*);
  void ExtractDoItData(G4ITStepProcessor*);

  void AddTrackID(G4Track*);

  void ResetLeadingTracks();

private:
  G4Scheduler();
  void Create();
  G4Scheduler(const G4Scheduler&);
  G4Scheduler& operator=(const G4Scheduler&);

  G4SchedulerMessenger* fSteppingMsg;

  static G4ThreadLocal G4Scheduler* fgScheduler;
  int fVerbose;
  bool fWhyDoYouStop;
  bool fInitialized;
  bool fRunning;

  int fNbTracks;
  int fNbSteps;
  int fMaxSteps;

  G4ITStepStatus fITStepStatus;

  // Time members
  double fTimeTolerance;
  double fGlobalTime;
  double fTmpGlobalTime;
  double fStartTime;
  double fEndTime;
  double fTmpEndTime;
  // fTmpEndTime : Stores the biggest end time here (for synchronizing tracks)
  double fPreviousTimeStep;
  int fZeroTimeCount;
  int fMaxNZeroTimeStepsAllowed;

  // Flags
  bool fComputeTimeStep;
  bool fComputeReaction;

  double fTimeStep; // The selected minimum time step

  // User steps
  bool fUsePreDefinedTimeSteps;
  double fDefaultMinTimeStep;
  std::map<double, double>* fpUserTimeSteps;
  // One can give time steps in respect to the global time
  mutable double fUserUpperTimeLimit;
  double fDefinedMinTimeStep;
  // selected user time step in respect to the global time
  bool fReachedUserTimeLimit; // if fMinTimeStep == the user time step

  double fTSTimeStep;
  // Time calculated by the time stepper in CalculateMinTimeStep()
  double fILTimeStep;
  // Time calculated by the interaction length methods
  // in ComputeInteractionLength()

  // std::map<G4Track*, G4TrackVectorHandle, compTrackPerID> fReactingTracks;
  G4ITReactionSet fReactionSet;
  std::vector<G4Track*> fLeadingTracks;

  bool fInteractionStep;
  // Flag : if the step is driven by the interaction with the matter and
  // NOT by the reaction between tracks

  G4ITTrackHolder& fTrackContainer;

  G4ITStepProcessor* fpStepProcessor;
  G4ITModelProcessor* fpModelProcessor;

  G4ITModelHandler* fpModelHandler;

  G4ITTrackingManager* fpTrackingManager;
  G4UserTimeStepAction* fpUserTimeStepAction;
  G4ITTrackingInteractivity* fpTrackingInteractivity;

  G4ITGun* fpGun;
  G4bool fContinue;
  G4bool fUseDefaultTimeSteps;
};

inline bool G4Scheduler::IsInitialized()
{
  return fInitialized;
}

inline G4ITModelHandler* G4Scheduler::GetModelHandler()
{
  return fpModelHandler;
}

inline void G4Scheduler::SetEndTime(const double __endtime)
{
  fEndTime = __endtime;
}

inline
void G4Scheduler::SetTimeSteps(std::map<double, double>* steps)
{
  fUsePreDefinedTimeSteps = true;
  fpUserTimeSteps = steps;
}

inline void G4Scheduler::AddTimeStep(double startingTime, double timeStep)
{
  if (fpUserTimeSteps == 0)
  {
    fpUserTimeSteps = new std::map<double, double>();
    fUsePreDefinedTimeSteps = true;
  }

  (*fpUserTimeSteps)[startingTime] = timeStep;
}

inline G4int G4Scheduler::GetNbSteps() const
{
  return fNbSteps;
}

inline void G4Scheduler::SetMaxNbSteps(G4int maxSteps)
{
  fMaxSteps = maxSteps;
}

inline G4int G4Scheduler::GetMaxNbSteps() const
{
  return fMaxSteps;
}

inline G4double G4Scheduler::GetStartTime() const
{
  return fStartTime;
}

inline G4double G4Scheduler::GetEndTime() const
{
  return fEndTime;
}

inline G4double G4Scheduler::GetTimeStep() const
{
  return fTimeStep;
}

inline void G4Scheduler::SetDefaultTimeStep(double timeStep)
{
  fDefaultMinTimeStep = timeStep;
}

inline G4double G4Scheduler::GetGlobalTime() const
{
  return fGlobalTime;
}

inline
void G4Scheduler::SetUserAction(G4UserTimeStepAction* userITAction)
{
  fpUserTimeStepAction = userITAction;
}

inline G4UserTimeStepAction* G4Scheduler::GetUserTimeStepAction() const
{
  return fpUserTimeStepAction;
}

inline void G4Scheduler::SetVerbose(int verbose)
{
  fVerbose = verbose;
}

inline int G4Scheduler::GetVerbose() const
{
  return fVerbose;
}

inline
void G4Scheduler::SetMaxZeroTimeAllowed(int maxTimeStepAllowed)
{
  fMaxNZeroTimeStepsAllowed = maxTimeStepAllowed;
}

inline int G4Scheduler::GetMaxZeroTimeAllowed() const
{
  return fMaxNZeroTimeStepsAllowed;
}

inline void G4Scheduler::SetTimeTolerance(double time)
{
  fTimeTolerance = time;
}

inline double G4Scheduler::GetTimeTolerance() const
{
  return fTimeTolerance;
}

inline G4double G4Scheduler::GetPreviousTimeStep() const
{
  return fPreviousTimeStep;
}

inline G4ITStepStatus G4Scheduler::GetStatus() const
{
  return fITStepStatus;
}

inline void G4Scheduler::Stop()
{
  fContinue = false;
}

inline G4ITTrackingInteractivity* G4Scheduler::GetInteractivity()
{
  return fpTrackingInteractivity;
}

inline void G4Scheduler::SetGun(G4ITGun* gun)
{
  fpGun = gun;
}

inline G4bool G4Scheduler::GetComputeTimeStepFlag() const
{
  return fComputeTimeStep;
}

inline void G4Scheduler::WhyDoYouStop()
{
  fWhyDoYouStop = true;
}

inline void G4Scheduler::UseDefaultTimeSteps(G4bool flag)
{
  fUseDefaultTimeSteps = flag;
}

inline G4bool G4Scheduler::AreDefaultTimeStepsUsed()
{
  return (fUseDefaultTimeSteps == false && fUsePreDefinedTimeSteps == false);
}

#endif
