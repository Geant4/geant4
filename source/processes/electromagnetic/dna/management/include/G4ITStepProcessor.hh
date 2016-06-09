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

#ifndef G4ITSTEPPROCESSOR_H
#define G4ITSTEPPROCESSOR_H

#include "G4ios.hh"                   // Include from 'system'
#include "globals.hh"                 // Include from 'global'
#include "Randomize.hh"               // Include from 'global'

#include "G4Navigator.hh"             // Include from 'geometry'
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

  
  G4Track* GetTrack()                                            {return fpTrack;}
  G4Step* GetStep()                                              {return fpStep;}
  const G4Step* GetStep() const                                  {return fpStep;}
  void SetStep(G4Step* val)                                      {fpStep = val;}

  inline G4TrackVector* GetSecondaries()                         {return fpSecondary;}
  inline void SetTrackingManager(G4ITTrackingManager* trackMan)  {fpTrackingManager = trackMan;}
  inline G4ITTrackingManager* GetTrackingManager()               {return fpTrackingManager;}


  virtual void Initialize(void* o = 0);
  void Initialize(G4Track*);
  
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

  void SetTrack(G4Track*);

  void GetProcessNumber();

  void SetupMembers();
  void ResetSecondaries();
  void InitDefineStep();
  void InitStepping(G4Track*);

  void SetInitialStep();

  void GetAtRestIL();
  void DoDefinePhysicalStepLength();
  G4StepStatus DoStepping();

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
  void SetNavigator(G4Navigator* value);
  G4double CalculateSafety();

  // Return the estimated safety value at the PostStepPoint
  void ApplyProductionCut(G4Track*);


  G4ITStepProcessor(const G4ITStepProcessor& other);
  G4ITStepProcessor& operator=(const G4ITStepProcessor& other);

private:

  G4ITTrackingManager* fpTrackingManager;

  // Member attributes
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

  G4ITTransportation* fpTransportation ;

  G4SelectedAtRestDoItVector* fpSelectedAtRestDoItVector;
  G4SelectedAlongStepDoItVector* fpSelectedAlongStepDoItVector;
  G4SelectedPostStepDoItVector* fpSelectedPostStepDoItVector;

  G4StepStatus* fpStepStatus;

  G4VITProcess* fpCurrentProcess;
  // The pointer to the process of which DoIt or
  // GetPhysicalInteractionLength has been just executed

  size_t MAXofAtRestLoops;
  size_t MAXofAlongStepLoops;
  size_t MAXofPostStepLoops;

  size_t fAtRestDoItProcTriggered;
  size_t fAlongStepDoItProcTriggered;
  size_t fPostStepDoItProcTriggered;

  G4int fN2ndariesAtRestDoIt;
  G4int fN2ndariesAlongStepDoIt;
  G4int fN2ndariesPostStepDoIt;
  // These are the numbers of secondaries generated by the process
  // just executed.

  //*********************************************************
  G4TrackVector*          fpSecondary ;

//  G4UserSteppingAction*   fpUserSteppingAction;

  G4bool                  KillVerbose;

  G4double*               fpPhysicalStep;    // Taken from G4TrackingInformation
  G4double                fPreviousStepSize;
  G4double                fTimeStep ;
  G4double                fPreviousStepTime ;

  G4VParticleChange*      fpParticleChange;
  G4Track*                fpTrack;
  G4IT*                   fpITrack ;
  G4TrackingInformation*  fpTrackingInfo ;
  G4Step*                 fpStep;
  G4StepPoint*            fpPreStepPoint;
  G4StepPoint*            fpPostStepPoint;

  G4VPhysicalVolume*      fpCurrentVolume;
  G4VSensitiveDetector*   fpSensitive;

  G4Navigator*            fpNavigator;
  G4int                   fStoreTrajectory;
  G4int                   verboseLevel;

  G4TouchableHandle* fpTouchableHandle;

  G4SteppingControl StepControlFlag;

  G4double kCarTolerance;
  // Cached geometrical tolerance on surface
  G4double proposedSafety;
  // This keeps the minimum safety value proposed by AlongStepGPILs.
  G4ThreeVector endpointSafOrigin;
  G4double endpointSafety;
  // To get the true safety value at the PostStepPoint, you have
  // to subtract the distance to 'endpointSafOrigin' from this value.
  G4double physIntLength;
  G4ForceCondition fCondition;
  G4GPILSelection  fGPILSelection;
  // Above three variables are for the method
  // DefinePhysicalStepLength(). To pass these information to
  // the method Verbose, they are kept at here. Need a more
  // elegant mechanism.


};

inline const G4Track* G4ITStepProcessor::GetTrack() const
{
  return fpTrack;
}

inline G4double G4ITStepProcessor::CalculateSafety()
{
  return std::max( endpointSafety -
		   (endpointSafOrigin - fpPostStepPoint->GetPosition()).mag(),
		   kCarTolerance );
}

inline void G4ITStepProcessor::SetNavigator(G4Navigator* value)
{
  fpNavigator = value;
}

inline void G4ITStepProcessor::CleanProcessor()
{
  fTimeStep = DBL_MAX ;
  fpPhysicalStep = 0 ;
  fpTrack = 0 ;
  fpITrack = 0;
  fpStep = 0 ;
  fpTrackingInfo = 0 ;

  fpParticleChange = 0;
  fpTrack = 0;
  fpITrack = 0;
  fpStep = 0;
  fpPreStepPoint = 0;
  fpPostStepPoint = 0;

  fpCurrentVolume = 0;
  fpSensitive = 0;

  fpSecondary = 0 ;
  fpAtRestDoItVector = 0;
  fpAlongStepDoItVector = 0;
  fpPostStepDoItVector= 0;

  fpAtRestGetPhysIntVector= 0;
  fpAlongStepGetPhysIntVector= 0;
  fpPostStepGetPhysIntVector= 0;

  fpTransportation = 0;

  fpSelectedAtRestDoItVector = 0;
  fpSelectedAlongStepDoItVector= 0;
  fpSelectedPostStepDoItVector= 0;

  fpStepStatus= 0;

  fpCurrentProcess= 0;
}

//______________________________________________________________________________
inline double G4ITStepProcessor::GetInteractionTime()
{
  return fTimeStep ;
}


#endif // G4ITSTEPPROCESSOR_H
