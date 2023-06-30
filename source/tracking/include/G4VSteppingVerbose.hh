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
// G4VSteppingVerbose
//
// Class description:
//
// This class manages the verbose outputs in G4SteppingManager.
// The instance should be a singleton. Users can inherit this
// class to make their own verbosity class.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------
#ifndef G4VSteppingVerbose_hh
#define G4VSteppingVerbose_hh 1

#include "G4ForceCondition.hh"  // enum 'track'
#include "G4GPILSelection.hh"  // enum 'track'
#include "G4StepStatus.hh"  // Include from 'track'
#include "G4TouchableHandle.hh"
#include "G4TrackVector.hh"  // Include from 'tracking'
#include "G4VProcess.hh"
#include "globals.hh"  // Include from 'global'

#include "trkgdefs.hh"

#include <vector>

class G4SteppingManager;
class G4Navigator;
class G4VPhysicalVolume;
class G4VSensitiveDetector;
class G4ProcessVector;
class G4SteppingManager;
class G4Track;
class G4UserSteppingAction;
class G4StepPoint;
class G4VParticleChange;

class G4VSteppingVerbose
{
 public:
  virtual ~G4VSteppingVerbose();

  // static methods to set/get the object's pointer

  static void SetInstance(G4VSteppingVerbose* Instance);
  static G4VSteppingVerbose* GetInstance();
  static G4VSteppingVerbose* GetMasterInstance();
  static G4int GetSilent();
  static void SetSilent(G4int fSilent);
  static G4int GetSilentStepInfo();
  static void SetSilentStepInfo(G4int fSilent);

  virtual G4VSteppingVerbose* Clone();

  // these method are invoked by G4SteppingManager

  virtual void NewStep() = 0;
  void CopyState();
  virtual void SetManager(G4SteppingManager* const);
  virtual void AtRestDoItInvoked() = 0;
  virtual void AlongStepDoItAllDone() = 0;
  virtual void PostStepDoItAllDone() = 0;
  virtual void AlongStepDoItOneByOne() = 0;
  virtual void PostStepDoItOneByOne() = 0;
  virtual void StepInfo() = 0;
  virtual void TrackingStarted() = 0;
  virtual void DPSLStarted() = 0;
  virtual void DPSLUserLimit() = 0;
  virtual void DPSLPostStep() = 0;
  virtual void DPSLAlongStep() = 0;
  virtual void VerboseTrack() = 0;
  virtual void VerboseParticleChange() = 0;

 protected:
  G4VSteppingVerbose();  // 'singleton'

  static G4ThreadLocal G4VSteppingVerbose* fInstance;  // pointer to the instance
  static G4VSteppingVerbose* fMasterInstance;  // pointer to the instance in master thread
  G4TRACKING_DLL static G4ThreadLocal G4int Silent;  // flag for verbosity
  G4TRACKING_DLL static G4ThreadLocal G4int SilentStepInfo;  // another flag for verbosity

  G4SteppingManager* fManager = nullptr;
  G4UserSteppingAction* fUserSteppingAction = nullptr;

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
  G4VProcess* fCurrentProcess = nullptr;  // The pointer to the process whose DoIt() or
                                          // GetPhysicalInteractionLength() has been just executed

  G4ProcessVector* fAtRestDoItVector = nullptr;
  G4ProcessVector* fAlongStepDoItVector = nullptr;
  G4ProcessVector* fPostStepDoItVector = nullptr;

  G4ProcessVector* fAtRestGetPhysIntVector = nullptr;
  G4ProcessVector* fAlongStepGetPhysIntVector = nullptr;
  G4ProcessVector* fPostStepGetPhysIntVector = nullptr;

  std::size_t MAXofAtRestLoops = 0;
  std::size_t MAXofAlongStepLoops = 0;
  std::size_t MAXofPostStepLoops = 0;

  G4double currentMinimumStep = 0.0;
  G4double numberOfInteractionLengthLeft = 0.0;

  std::size_t fAtRestDoItProcTriggered = 0;
  std::size_t fAlongStepDoItProcTriggered = 0;
  std::size_t fPostStepDoItProcTriggered = 0;

  G4int fN2ndariesAtRestDoIt = 0;
  G4int fN2ndariesAlongStepDoIt = 0;
  G4int fN2ndariesPostStepDoIt = 0;
  // These are the numbers of secondaries generated by the process
  // just executed

  G4Navigator* fNavigator = nullptr;

  G4int verboseLevel = 0;

  using G4SelectedAtRestDoItVector = std::vector<G4int>;
  using G4SelectedAlongStepDoItVector = std::vector<G4int>;
  using G4SelectedPostStepDoItVector = std::vector<G4int>;

  G4SelectedAtRestDoItVector* fSelectedAtRestDoItVector = nullptr;
  G4SelectedAlongStepDoItVector* fSelectedAlongStepDoItVector = nullptr;
  G4SelectedPostStepDoItVector* fSelectedPostStepDoItVector = nullptr;

  G4double fPreviousStepSize = 0.0;

  G4TouchableHandle fTouchableHandle;

  G4SteppingControl StepControlFlag = NormalCondition;

  G4double physIntLength = 0.0;
  G4ForceCondition fCondition = InActivated;
  G4GPILSelection fGPILSelection = NotCandidateForSelection;
  // Above three variables are for the method DefinePhysicalStepLength().
  // To pass this information to the method Verbose(), they are kept at
  // here. Need a more elegant mechanism
};

#endif
