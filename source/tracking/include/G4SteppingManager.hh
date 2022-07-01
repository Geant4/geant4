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
// G4SteppingManager
//
// Class description:
//
// This is the class which plays an essential role in tracking particles.
// It takes cares of all message passing between objects in the different
// categories (for example, geometry - including transportation, interactions
// in matter, etc). Its public method 'stepping' steers to step the particle.
// Only used within the Geant4 kernel.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// History:
//   12.03.1998, H.Kurashige - modified for use of G4ParticleChange
// --------------------------------------------------------------------
#ifndef G4SteppingManager_hh
#define G4SteppingManager_hh 1

#include <iomanip>                     // Include from 'system'
#include <vector>                      // Include from 'system'
#include "globals.hh"                  // Include from 'global'
#include "Randomize.hh"                // Include from 'global'

#include "G4Navigator.hh"              // Include from 'geometry'
#include "G4LogicalVolume.hh"          // Include from 'geometry'
#include "G4VPhysicalVolume.hh"        // Include from 'geometry'
#include "G4ProcessManager.hh"         // Include from 'processes'
#include "G4NoProcess.hh"              // Include from 'processes'

#include "G4Track.hh"                  // Include from 'tracking'
#include "G4TrackVector.hh"            // Include from 'tracking'
#include "G4TrackStatus.hh"            // Include from 'tracking'
#include "G4StepStatus.hh"             // Include from 'tracking'
#include "G4UserSteppingAction.hh"     // Include from 'tracking'
#include "G4Step.hh"                   // Include from 'tracking'
#include "G4StepPoint.hh"              // Include from 'tracking'
#include "G4VSteppingVerbose.hh"       // Include from 'tracking'
#include "G4TouchableHandle.hh"        // Include from 'geometry'
#include "G4TouchableHistoryHandle.hh" // Include from 'geometry'

using G4SelectedAtRestDoItVector = std::vector<G4int>;
using G4SelectedAlongStepDoItVector = std::vector<G4int>;
using G4SelectedPostStepDoItVector = std::vector<G4int>;

class G4VSensitiveDetector;

class G4SteppingManager 
{
 public:
  using ProfilerConfig = G4Step::ProfilerConfig;

 public:
    // Constructor/Destructor

    G4SteppingManager();
      // SteppingManger should be dynamically allocated, therefore 
      // you need to invoke new() when you call this constructor.
      // "Secodary track vector" will be dynamically created by this 
      // constructor. G4UserSteppingAction will be also created 
      // in this constructor, and "this" pointer will be passed to 
      // G4UserSteppingAction.

   ~G4SteppingManager();

    // Get/Set functions

    const G4TrackVector* GetSecondary() const;
    void SetUserAction(G4UserSteppingAction* apAction);
    G4Track* GetTrack() const;
    void SetVerboseLevel(G4int vLevel);
    void SetVerbose(G4VSteppingVerbose*);
    G4Step* GetStep() const;
    void SetNavigator(G4Navigator* value);

    // Other member functions

    G4StepStatus Stepping();
      // Steers to move the give particle from the TrackingManger by one Step.

    void SetInitialStep(G4Track* valueTrack);
      // Sets up initial track information (enegry, position, etc) to 
      // the PreStepPoint of the G4Step. This funciton has to be called 
      // just once before the stepping loop in the "TrackingManager".

    void GetProcessNumber();

    // Get methods

    G4double GetPhysicalStep();
    G4double GetGeometricalStep();
    G4double GetCorrectedStep();
    G4bool GetPreStepPointIsGeom();
    G4bool GetFirstStep();
    G4StepStatus GetfStepStatus();
    G4double GetTempInitVelocity();
    G4double GetTempVelocity();
    G4double GetMass();
    G4double GetsumEnergyChange();
    G4VParticleChange* GetfParticleChange();
    G4Track* GetfTrack();
    G4TrackVector* GetfSecondary();
    G4Step* GetfStep();
    G4StepPoint* GetfPreStepPoint();
    G4StepPoint* GetfPostStepPoint();
    G4VPhysicalVolume* GetfCurrentVolume();
    G4VSensitiveDetector* GetfSensitive();
    G4VProcess* GetfCurrentProcess();
    G4ProcessVector* GetfAtRestDoItVector();
    G4ProcessVector* GetfAlongStepDoItVector();
    G4ProcessVector* GetfPostStepDoItVector();
    G4ProcessVector* GetfAlongStepGetPhysIntVector();
    G4ProcessVector* GetfPostStepGetPhysIntVector();
    G4ProcessVector* GetfAtRestGetPhysIntVector();
    G4double GetcurrentMinimumStep();
    G4double GetnumberOfInteractionLengthLeft();
    std::size_t GetfAtRestDoItProcTriggered();
    std::size_t GetfAlongStepDoItProcTriggered();
    std::size_t GetfPostStepDoItProcTriggered();
    G4int GetfN2ndariesAtRestDoIt();
    G4int GetfN2ndariesAlongStepDoIt();
    G4int GetfN2ndariesPostStepDoIt();
    G4Navigator* GetfNavigator();
    G4int GetverboseLevel();
    std::size_t GetMAXofAtRestLoops();
    std::size_t GetMAXofAlongStepLoops();
    std::size_t GetMAXofPostStepLoops();
    G4SelectedAtRestDoItVector* GetfSelectedAtRestDoItVector();
    G4SelectedAlongStepDoItVector* GetfSelectedAlongStepDoItVector();
    G4SelectedPostStepDoItVector* GetfSelectedPostStepDoItVector();
    G4double GetfPreviousStepSize();
    const G4TouchableHandle& GetTouchableHandle();
    G4SteppingControl GetStepControlFlag();
    G4UserSteppingAction* GetUserAction();
    G4double GetphysIntLength();
    G4ForceCondition GetfCondition();
    G4GPILSelection GetfGPILSelection();

  private:

    // Member functions

    void DefinePhysicalStepLength();
      // Calculate corresponding physical length from the mean free path 
      // left for each discrete physics process. The minimum allowable
      // step for each continuous process will be also calculated.
    void InvokeAtRestDoItProcs();
    void InvokeAlongStepDoItProcs();
    void InvokePostStepDoItProcs();
    void InvokePSDIP(size_t); // 
    G4int ProcessSecondariesFromParticleChange();
    G4double CalculateSafety();
      // Return the estimated safety value at the PostStepPoint

    // Member data 

    static const size_t SizeOfSelectedDoItVector = 100;

    G4bool KillVerbose = false;

    G4UserSteppingAction* fUserSteppingAction = nullptr;

    G4VSteppingVerbose* fVerbose = nullptr;

    G4double PhysicalStep = 0.0;
    G4double GeometricalStep = 0.0;
    G4double CorrectedStep = 0.0;
    G4bool PreStepPointIsGeom = false;
    G4bool FirstStep = false;
    G4StepStatus fStepStatus = fUndefined;

    G4double TempInitVelocity = 0.0;
    G4double TempVelocity = 0.0;
    G4double Mass = 0.0;

    G4double sumEnergyChange = 0.0;

    G4VParticleChange* fParticleChange = nullptr;
    G4Track* fTrack = nullptr;
    G4TrackVector* fSecondary = nullptr;
    G4Step* fStep = nullptr;
    G4StepPoint* fPreStepPoint = nullptr;
    G4StepPoint* fPostStepPoint = nullptr;

    G4VPhysicalVolume* fCurrentVolume = nullptr;
    G4VSensitiveDetector* fSensitive = nullptr;
    G4VProcess* fCurrentProcess = nullptr;
      // The pointer to the process of which DoIt() or
      // GetPhysicalInteractionLength() has been just executed.


    G4ProcessVector* fAtRestDoItVector = nullptr;
    G4ProcessVector* fAlongStepDoItVector = nullptr;
    G4ProcessVector* fPostStepDoItVector = nullptr;

    G4ProcessVector* fAtRestGetPhysIntVector = nullptr;
    G4ProcessVector* fAlongStepGetPhysIntVector = nullptr;
    G4ProcessVector* fPostStepGetPhysIntVector = nullptr;

    std::size_t MAXofAtRestLoops = 0;
    std::size_t MAXofAlongStepLoops = 0;
    std::size_t MAXofPostStepLoops = 0;

    std::size_t fAtRestDoItProcTriggered = 0;
    std::size_t fAlongStepDoItProcTriggered = 0;
    std::size_t fPostStepDoItProcTriggered = 0;

    G4int fN2ndariesAtRestDoIt = 0;
    G4int fN2ndariesAlongStepDoIt = 0;
    G4int fN2ndariesPostStepDoIt = 0;
      // These are the numbers of secondaries generated by the process
      // just executed.

    G4Navigator* fNavigator = nullptr;

    G4int verboseLevel = 0;

    G4SelectedAtRestDoItVector* fSelectedAtRestDoItVector = nullptr;
    G4SelectedAlongStepDoItVector* fSelectedAlongStepDoItVector = nullptr;
    G4SelectedPostStepDoItVector* fSelectedPostStepDoItVector = nullptr;

    G4double fPreviousStepSize = 0.0;

    G4TouchableHandle fTouchableHandle;

    G4SteppingControl StepControlFlag = NormalCondition;

    G4double kCarTolerance = 0.0;
      // Cached geometrical tolerance on surface
    G4double proposedSafety = 0.0;
      // This keeps the minimum safety value proposed by AlongStepGPILs.
    G4ThreeVector endpointSafOrigin;
    G4double endpointSafety = 0.0;
      // To get the true safety value at the PostStepPoint, you have
      // to subtract the distance to 'endpointSafOrigin' from this value.
    G4double physIntLength = 0.0;
    G4ForceCondition fCondition = InActivated;
    G4GPILSelection  fGPILSelection = NotCandidateForSelection;
      // Above three variables are for the method 
      // DefinePhysicalStepLength(). To pass these information to
      // the method Verbose, they are kept at here. Need a more 
      // elegant mechanism.

    G4NoProcess const* fNoProcess = nullptr;
      // Used in the InvokeAtRestDoItProcs() method to flag the process
      // of any stable ion at rest. 
};

//*******************************************************************
//
// Inline functions
//
//*******************************************************************

  inline G4double G4SteppingManager::GetPhysicalStep()
  {
    return PhysicalStep;
  }

  inline  G4double G4SteppingManager::GetGeometricalStep()
  {
    return GeometricalStep;
  }

  inline  G4double G4SteppingManager::GetCorrectedStep()
  {
    return CorrectedStep;
  }
  
  inline G4bool G4SteppingManager::GetPreStepPointIsGeom()
  {
    return PreStepPointIsGeom;
  }

  inline G4bool G4SteppingManager::GetFirstStep()
  {
    return FirstStep;
  }

  inline G4StepStatus G4SteppingManager::GetfStepStatus()
  {
    return fStepStatus;
  }

  inline G4double G4SteppingManager::GetTempInitVelocity()
  {
    return TempInitVelocity;
  }

  inline G4double G4SteppingManager::GetTempVelocity()
  {
    return TempVelocity;
  }

  inline G4double G4SteppingManager::GetMass()
  {
    return Mass;
  }

  inline G4double G4SteppingManager::GetsumEnergyChange()
  {
    return sumEnergyChange;
  }

  inline G4VParticleChange* G4SteppingManager::GetfParticleChange()
  {
    return fParticleChange;
  }

  inline G4Track* G4SteppingManager::GetfTrack()
  {
    return fTrack;
  }

  inline G4TrackVector* G4SteppingManager::GetfSecondary()
  {
    return fStep->GetfSecondary();
  }

  inline G4Step* G4SteppingManager::GetfStep()
  {
    return fStep;
  }

  inline G4StepPoint* G4SteppingManager::GetfPreStepPoint()
  {
    return fPreStepPoint;
  }

  inline G4StepPoint* G4SteppingManager::GetfPostStepPoint()
  {
    return fPostStepPoint;
  }

  inline G4VPhysicalVolume* G4SteppingManager::GetfCurrentVolume()
  {
    return fCurrentVolume;
  }

  inline G4VSensitiveDetector* G4SteppingManager::GetfSensitive()
  {
    return fSensitive;
  }

  inline G4VProcess* G4SteppingManager::GetfCurrentProcess()
  {
    return fCurrentProcess;
  }

  inline G4ProcessVector* G4SteppingManager::GetfAtRestDoItVector()
  {
    return fAtRestDoItVector;
  }

  inline G4ProcessVector* G4SteppingManager::GetfAlongStepDoItVector()
  {
    return fAlongStepDoItVector;
  }

  inline G4ProcessVector* G4SteppingManager::GetfPostStepDoItVector()
  {
    return fPostStepDoItVector;
  }

  inline G4ProcessVector* G4SteppingManager::GetfAtRestGetPhysIntVector()
  {
    return fAtRestGetPhysIntVector;
  }

  inline G4ProcessVector* G4SteppingManager::GetfAlongStepGetPhysIntVector()
  {
    return fAlongStepGetPhysIntVector;
  }

  inline G4ProcessVector* G4SteppingManager::GetfPostStepGetPhysIntVector()
  {
    return fPostStepGetPhysIntVector;
  }

  inline size_t G4SteppingManager::GetMAXofAtRestLoops()
  {
    return MAXofAtRestLoops;
  }

  inline size_t G4SteppingManager::GetMAXofAlongStepLoops()
  {
    return MAXofAlongStepLoops;
  }

  inline size_t G4SteppingManager::GetMAXofPostStepLoops()
  {
    return MAXofPostStepLoops;
  }

  inline size_t G4SteppingManager::GetfAtRestDoItProcTriggered()
  {
    return fAtRestDoItProcTriggered;
  }

  inline size_t G4SteppingManager::GetfAlongStepDoItProcTriggered()
  {
    return fAtRestDoItProcTriggered;
  }

  inline size_t G4SteppingManager::GetfPostStepDoItProcTriggered()
  {
    return fPostStepDoItProcTriggered;
  }

  inline G4int G4SteppingManager::GetfN2ndariesAtRestDoIt()
  {
    return fN2ndariesAtRestDoIt;
  }

  inline G4int G4SteppingManager::GetfN2ndariesAlongStepDoIt()
  {
    return fN2ndariesAlongStepDoIt;
  }

  inline G4int G4SteppingManager::GetfN2ndariesPostStepDoIt()
  {
    return fN2ndariesPostStepDoIt;
  }

  inline G4Navigator* G4SteppingManager::GetfNavigator()
  {
    return fNavigator;
  }

  inline G4int G4SteppingManager::GetverboseLevel()
  {
    return verboseLevel;
  }

  inline G4SelectedAtRestDoItVector*
  G4SteppingManager::GetfSelectedAtRestDoItVector()
  {
    return fSelectedAtRestDoItVector;
  }

  inline G4SelectedAlongStepDoItVector*
  G4SteppingManager::GetfSelectedAlongStepDoItVector()
  {
    return fSelectedAlongStepDoItVector;
  }

  inline G4SelectedPostStepDoItVector*
  G4SteppingManager::GetfSelectedPostStepDoItVector()
  {
    return fSelectedPostStepDoItVector;
  }

  inline G4double G4SteppingManager::GetfPreviousStepSize()
  {
    return fPreviousStepSize;
  }

  inline const G4TouchableHandle& G4SteppingManager::GetTouchableHandle()
  {
    return fTouchableHandle;
  }

  inline G4SteppingControl G4SteppingManager::GetStepControlFlag()
  {
    return StepControlFlag;
  }

  inline G4double G4SteppingManager::GetphysIntLength()
  {
    return physIntLength;
  }

  inline G4ForceCondition G4SteppingManager::GetfCondition()
  {
    return fCondition;
  }

  inline G4GPILSelection  G4SteppingManager::GetfGPILSelection()
  {
    return fGPILSelection;
  }

  inline const G4TrackVector* G4SteppingManager::GetSecondary() const
  {
    return fStep->GetSecondary(); 
  }

  inline void G4SteppingManager::SetNavigator(G4Navigator* value)
  {
    fNavigator = value; 
  }

  inline void G4SteppingManager::SetUserAction(G4UserSteppingAction* apAction)
  {
    fUserSteppingAction = apAction;
  }

  inline G4UserSteppingAction* G4SteppingManager::GetUserAction()
  {
    return fUserSteppingAction;
  }

  inline G4Track* G4SteppingManager::GetTrack() const
  {
    return fTrack;
  }

  inline void G4SteppingManager::SetVerboseLevel(G4int vLevel)
  {
    verboseLevel = vLevel; 
  }

  inline void G4SteppingManager::SetVerbose(G4VSteppingVerbose* yourVerbose)
  {
    fVerbose = yourVerbose;
  }

  inline G4Step* G4SteppingManager::GetStep() const
  {
    return fStep;
  }

  inline G4double G4SteppingManager::CalculateSafety()
  {
    return std::max( endpointSafety -
                (endpointSafOrigin - fPostStepPoint->GetPosition()).mag(),
                kCarTolerance );
  }

#endif
