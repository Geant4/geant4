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
#include "G4VScavengerMaterial.hh"

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
    G4bool operator()(G4Track* rhs, G4Track* lhs) const
    {
      return rhs->GetTrackID() < lhs->GetTrackID();
    }
  };
#endif

/**
 * G4Scheduler synchronizes (in time) track stepping
 */
class G4Scheduler :
    public G4VScheduler,
    public G4VStateDependent
{
protected:
  ~G4Scheduler() override;

public:
  G4Scheduler(const G4Scheduler&) = delete;
  G4Scheduler& operator=(const G4Scheduler&) = delete;

  static G4Scheduler* Instance();
  /** DeleteInstance should be used instead
   * of the destructor
   */
  static void DeleteInstance();
  G4bool Notify(G4ApplicationState requestedState) override;

  void RegisterModel(G4VITStepModel*, G4double) override;

  void Initialize() override;
  void ForceReinitialization();
  inline G4bool IsInitialized();
  inline G4bool IsRunning() override{return fRunning;}
  void Reset() override;
  void Process() override;
  void ClearList();

  inline void SetGun(G4ITGun*) override;
  inline G4ITGun* GetGun();

  inline void Stop();
  void Clear();

  // To be called only in UserReactionAction::EndProcessing()
  // after fRunning flag has been turned off.
  // This is not done automatically before UserReactionAction::EndProcessing()
  // is called in case one would like to access some track information
  void EndTracking();

  void SetEndTime(const G4double) override;

  /* Two tracks below the time tolerance are supposed to be
   * in the same time slice
   */
  inline void SetTimeTolerance(G4double) override;
  inline G4double GetTimeTolerance() const override;

  inline void SetMaxZeroTimeAllowed(G4int) override;
  inline G4int GetMaxZeroTimeAllowed() const override;

  inline G4ITModelHandler* GetModelHandler() override;

  inline void SetTimeSteps(std::map<G4double, G4double>*) override;
  inline void AddTimeStep(G4double, G4double) override;
  inline void SetDefaultTimeStep(G4double) override;
  G4double GetLimitingTimeStep() const override;
  inline G4int GetNbSteps() const override;
  inline void SetMaxNbSteps(G4int) override;
  inline G4int GetMaxNbSteps() const override;
  inline G4double GetStartTime() const override;
  inline G4double GetEndTime() const override;
  inline G4double GetTimeStep() const override;
  inline G4double GetPreviousTimeStep() const override;
  inline G4double GetGlobalTime() const override;
  inline void SetUserAction(G4UserTimeStepAction*) override;
  inline G4UserTimeStepAction* GetUserTimeStepAction() const override;

  // To use with transportation only, no reactions
  inline void UseDefaultTimeSteps(G4bool);
  inline G4bool AreDefaultTimeStepsUsed();

  inline G4ITStepStatus GetStatus() const;

  /* 1 : Reaction information
   * 2 : (1) + time step information
   * 3 : (2) + step info for individual tracks
   * 4 : (2) + trackList processing info + pushed and killed track info
   */
  inline void SetVerbose(G4int) override;
  
  inline G4int GetVerbose() const;

  inline void WhyDoYouStop();

  void SetInteractivity(G4ITTrackingInteractivity*) override;
  inline G4ITTrackingInteractivity* GetInteractivity() override;

  virtual size_t GetNTracks();

  void GetCollisionType(G4String& interactionType);

  void AddWatchedTime(G4double time)
  {
    fWatchedTimes.insert(time);
  }

  G4double GetNextWatchedTime() const;

  inline void SetMaxTimeStep(G4double maxTimeStep)
  {
    fMaxTimeStep = maxTimeStep;
  }

  inline G4double GetMaxTimeStep() const
  {
    return fMaxTimeStep;
  }

  inline G4VScavengerMaterial* GetScavengerMaterial() const
  {
      return fpUserScavenger.get();
  }
  inline void SetScavengerMaterial(std::unique_ptr<G4VScavengerMaterial> scavengerMaterial)
  {
      fpUserScavenger = std::move(scavengerMaterial);
  }

protected:

  void DoProcess();
  void SynchronizeTracks();
  void Stepping();

  void FindUserPreDefinedTimeStep();

  G4bool CanICarryOn();

  void PrintWhyDoYouStop();

private:
  G4Scheduler();
  void Create();

  G4SchedulerMessenger* fpMessenger;

  static G4ThreadLocal G4Scheduler* fgScheduler;
  G4int fVerbose;
  G4bool fWhyDoYouStop;
  G4bool fInitialized;
  G4bool fRunning;
  G4bool fContinue;

  G4int fNbSteps;
  G4int fMaxSteps;

  G4ITStepStatus fITStepStatus;

  // Time members
  G4bool fUseDefaultTimeSteps;
  G4double fTimeTolerance;
  G4double fGlobalTime;
  G4double fTmpGlobalTime;
  G4double fStartTime;
  G4double fStopTime;
  G4double fEndTime;
  G4double fPreviousTimeStep;
  G4int fZeroTimeCount;
  G4int fMaxNZeroTimeStepsAllowed;

  G4double fTimeStep; // The selected minimum time step
  G4double fMaxTimeStep;

  // User steps
  G4bool fUsePreDefinedTimeSteps;
  G4double fDefaultMinTimeStep;
  std::map<G4double, G4double>* fpUserTimeSteps;
  // One can give time steps in respect to the global time
  mutable G4double fUserUpperTimeLimit;
  G4double fDefinedMinTimeStep;
  // selected user time step in respect to the global time
  G4bool fReachedUserTimeLimit; // if fMinTimeStep == the user time step

  std::set<G4double> fWatchedTimes;

  G4UserTimeStepAction* fpUserTimeStepAction;

  std::unique_ptr<G4VScavengerMaterial> fpUserScavenger;

  // ==========================================
  // TO BE REMOVED
  G4ITStepProcessor* fpStepProcessor;
  G4ITModelProcessor* fpModelProcessor;
  G4ITTrackingManager* fpTrackingManager;
  G4ITTrackingInteractivity* fpTrackingInteractivity;
  G4ITReactionSet* fReactionSet;
  G4ITTrackHolder& fTrackContainer;
  G4ITModelHandler* fpModelHandler;
  // ==========================================

  G4double fTSTimeStep;
  // Time calculated by the time stepper in CalculateMinTimeStep()
  G4double fILTimeStep;
  // Time calculated by the interaction length methods
  // in ComputeInteractionLength()

  G4bool fInteractionStep;
  // Flag : if the step is driven by the interaction with the matter and
  // NOT by the reaction between tracks

  G4ITGun* fpGun;

  // ==========================================
  //Hoang
  bool fResetScavenger;
public:
  void ResetScavenger(bool);
};

inline G4bool G4Scheduler::IsInitialized()
{
  return fInitialized;
}

inline G4ITModelHandler* G4Scheduler::GetModelHandler()
{
  return fpModelHandler;
}

inline void G4Scheduler::SetEndTime(const G4double __endtime)
{
  fEndTime = __endtime;
}

inline
void G4Scheduler::SetTimeSteps(std::map<G4double, G4double>* steps)
{
  fUsePreDefinedTimeSteps = true;
  fpUserTimeSteps = steps;
}

inline void G4Scheduler::AddTimeStep(G4double startingTime, G4double timeStep)
{
  if (fpUserTimeSteps == nullptr)
  {
    fpUserTimeSteps = new std::map<G4double, G4double>();
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

inline void G4Scheduler::SetDefaultTimeStep(G4double timeStep)
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

inline void G4Scheduler::SetVerbose(G4int verbose)
{
  fVerbose = verbose;
}

inline G4int G4Scheduler::GetVerbose() const
{
  return fVerbose;
}

inline
void G4Scheduler::SetMaxZeroTimeAllowed(G4int maxTimeStepAllowed)
{
  fMaxNZeroTimeStepsAllowed = maxTimeStepAllowed;
}

inline G4int G4Scheduler::GetMaxZeroTimeAllowed() const
{
  return fMaxNZeroTimeStepsAllowed;
}

inline void G4Scheduler::SetTimeTolerance(G4double time)
{
  fTimeTolerance = time;
}

inline G4double G4Scheduler::GetTimeTolerance() const
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

inline G4ITGun* G4Scheduler::GetGun()
{
  return fpGun;
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
  return (!fUseDefaultTimeSteps && !fUsePreDefinedTimeSteps);
}

inline void G4Scheduler::ResetScavenger(bool value)
{
    fResetScavenger = value;
}

#endif
