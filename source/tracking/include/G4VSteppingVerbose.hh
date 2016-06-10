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
// $Id: G4VSteppingVerbose.hh 66872 2013-01-15 01:25:57Z japost $
//
//---------------------------------------------------------------
//
// G4VSteppingVerbose.hh
//
// class description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   The instance should be singleton. Users can inherit this 
//   class to make their own verbosing class.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4VSteppingVerbose_h
#define G4VSteppingVerbose_h

class G4SteppingManager;

#include "globals.hh"                 // Include from 'global'
#include <vector>

class G4Navigator;
class G4VPhysicalVolume;
class G4VSensitiveDetector;
#include "G4VProcess.hh"
class G4ProcessVector;
class G4SteppingManager;              // Include from 'tracking'
class G4Track;
#include "G4TrackVector.hh"           // Include from 'tracking'
#include "G4StepStatus.hh"            // Include from 'track'
class G4UserSteppingAction;
class G4StepPoint;
#include "G4TouchableHandle.hh"
class G4VParticleChange;
#include "G4ForceCondition.hh"  //enum 'track'
#include "G4GPILSelection.hh"   //enum 'track'


class G4VSteppingVerbose
{
// Constructor/Destructor
protected:    // to force 'singleton'
  G4VSteppingVerbose();  
public:  
  virtual ~G4VSteppingVerbose();
  //
protected:  
  static G4ThreadLocal G4VSteppingVerbose* fInstance;// pointer to the instance 
  static G4ThreadLocal G4int Silent; //flag for verbosity
  static G4ThreadLocal G4int SilentStepInfo; //another flag for verbosity
public:   // with description
// static methods to set/get the object's pointer 
  static void SetInstance(G4VSteppingVerbose* Instance);
  static G4VSteppingVerbose* GetInstance();
  static G4int GetSilent();
  static void SetSilent(G4int fSilent);
  static G4int GetSilentStepInfo();
  static void SetSilentStepInfo(G4int fSilent);
// these method are invoked in the SteppingManager 
  virtual void NewStep() = 0;
  void CopyState();
  void SetManager(G4SteppingManager* const);
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
  // Member data

protected:
  G4SteppingManager* fManager;
  
  G4UserSteppingAction* fUserSteppingAction;
  
  G4double PhysicalStep;
  G4double GeometricalStep;
  G4double CorrectedStep;
  G4bool PreStepPointIsGeom;
  G4bool FirstStep;
  G4StepStatus fStepStatus;

  G4double TempInitVelocity;
  G4double TempVelocity;
  G4double Mass;

  G4double sumEnergyChange;

  G4VParticleChange* fParticleChange;
  G4Track* fTrack;
  G4TrackVector* fSecondary;
  G4Step* fStep;
  G4StepPoint* fPreStepPoint;
  G4StepPoint* fPostStepPoint;

  G4VPhysicalVolume* fCurrentVolume;
  G4VSensitiveDetector* fSensitive;
  G4VProcess* fCurrentProcess;
  // The pointer to the process of which DoIt or
  // GetPhysicalInteractionLength has been just executed.


  G4ProcessVector* fAtRestDoItVector;
  G4ProcessVector* fAlongStepDoItVector;
  G4ProcessVector* fPostStepDoItVector;

  G4ProcessVector* fAtRestGetPhysIntVector;
  G4ProcessVector* fAlongStepGetPhysIntVector;
  G4ProcessVector* fPostStepGetPhysIntVector;

  size_t MAXofAtRestLoops;
  size_t MAXofAlongStepLoops;
  size_t MAXofPostStepLoops;
  
  G4double currentMinimumStep;
  G4double numberOfInteractionLengthLeft;

  size_t fAtRestDoItProcTriggered;
  size_t fAlongStepDoItProcTriggered;
  size_t fPostStepDoItProcTriggered;

  G4int fN2ndariesAtRestDoIt;
  G4int fN2ndariesAlongStepDoIt;
  G4int fN2ndariesPostStepDoIt;
      // These are the numbers of secondaries generated by the process
      // just executed.

  G4Navigator *fNavigator;

  G4int verboseLevel;

  typedef std::vector<G4int> 
             G4SelectedAtRestDoItVector;
  typedef std::vector<G4int> 
             G4SelectedAlongStepDoItVector;
  typedef std::vector<G4int>
             G4SelectedPostStepDoItVector;
  G4SelectedAtRestDoItVector* fSelectedAtRestDoItVector;
  G4SelectedAlongStepDoItVector* fSelectedAlongStepDoItVector;
  G4SelectedPostStepDoItVector* fSelectedPostStepDoItVector;

  G4double   fPreviousStepSize;

  G4TouchableHandle fTouchableHandle;

  G4SteppingControl StepControlFlag;

  G4double physIntLength;
  G4ForceCondition fCondition;
  G4GPILSelection  fGPILSelection;
      // Above three variables are for the method 
      // DefinePhysicalStepLength(). To pass these information to
      // the method Verbose, they are kept at here. Need a more 
      // elegant mechanism.


};
#endif
