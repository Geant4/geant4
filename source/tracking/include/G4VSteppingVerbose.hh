//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VSteppingVerbose.hh,v 1.15 2003-06-16 17:13:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
#include "G4TrackVector.hh"            // Include from 'tracking'
#include "G4StepStatus.hh"            // Include from 'track'
class G4UserSteppingAction;
class G4StepPoint;
#include "G4TouchableHandle.hh"
class G4VParticleChange;
#include "G4ForceCondition.hh"  //enum 'track'
#include "G4GPILSelection.hh"   //enum 'track'


class G4VSteppingVerbose{

// Constructor/Destructor
protected:    // to force 'singleton'
  G4VSteppingVerbose();  
public:  
  virtual ~G4VSteppingVerbose();
  //
public:   // with description
// static methods to set/get the object's pointer 
  static void SetInstance(G4VSteppingVerbose* Instance);
  static G4VSteppingVerbose* GetInstance();
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

protected:  // pointer to the instance 
  static G4VSteppingVerbose* fInstance;
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
