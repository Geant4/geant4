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
/*
 * G4VITSteppingVerbose.hh
 *
 *  Created on: Jun 22, 2015
 *      Author: mkaramit
 */

#ifndef SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4VITSTEPPINGVERBOSE_HH_
#define SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4VITSTEPPINGVERBOSE_HH_

#include "globals.hh"                 // Include from 'global'
#include <vector>

class G4Navigator;
class G4VPhysicalVolume;
class G4VSensitiveDetector;
#include "G4VProcess.hh"
class G4ProcessVector;
class G4Track;

//#include "G4TrackVector.hh"           // Include from 'tracking'
#include "G4StepStatus.hh"            // Include from 'track'
class G4UserSteppingAction;
class G4StepPoint;
#include "G4TouchableHandle.hh"

#include "G4ForceCondition.hh"  //enum 'track'
#include "G4GPILSelection.hh"   //enum 'track'
#include "G4ITStepProcessor.hh"
#include "G4VITProcess.hh"

class G4VParticleChange;
class G4ITStepProcessorState;
//class ProcessGeneralInfo;
//class G4VPhysicalVolume;
//class G4ProcessVector;

#include <G4UImessenger.hh>

class G4UIcmdWithAnInteger;

class G4VITSteppingVerbose : G4UImessenger
{
public:
  G4VITSteppingVerbose();
  virtual ~G4VITSteppingVerbose();

public:

  virtual void TrackingStarted(G4Track* track);
  virtual void TrackingEnded(G4Track* track);

  virtual void DoItStarted() = 0;
  virtual void PreStepVerbose(G4Track*) = 0;
  virtual void PostStepVerbose(G4Track*) = 0;

  // these methods are invoked in the SteppingManager
  virtual void NewStep() = 0;
  void CopyState();

  virtual void StepInfoForLeadingTrack() = 0;

  virtual void AtRestDoItInvoked() = 0;
  virtual void AtRestDoItOneByOne() = 0;

  virtual void PostStepDoItAllDone() = 0;
  virtual void PostStepDoItOneByOne() = 0;

  virtual void AlongStepDoItAllDone() = 0;
  virtual void AlongStepDoItOneByOne() = 0;

  virtual void StepInfo() = 0;
  virtual void DPSLStarted() = 0;
  virtual void DPSLUserLimit() = 0;
  virtual void DPSLPostStep() = 0;
  virtual void DPSLAlongStep() = 0;
  virtual void VerboseTrack() = 0;
  virtual void VerboseParticleChange() = 0;

  //____________________________________________________________________________

  inline void SetVerbose(int flag)
  {
    fVerboseLevel = flag;
  }

  inline G4int GetVerbose()
  {
    return fVerboseLevel;
  }

  //____________________________________________________________________________

  virtual void SetNewValue(G4UIcommand * command,
                           G4String newValue);

  virtual G4String GetCurrentValue(G4UIcommand * command);


  //____________________________________________________________________________

  void SetStepProcessor(const G4ITStepProcessor* stepProcessor)
  {
    this->fpStepProcessor = stepProcessor;
  }

  void TrackBanner(G4Track* track, const G4String& message);

protected:
  const G4ITStepProcessor* fpStepProcessor;

  G4UIcmdWithAnInteger* fpVerboseUI;
  G4ITStepProcessorState* fpState;
  const ProcessGeneralInfo* fpProcessGeneralInfo;

  G4double PhysicalStep;
  G4StepStatus fStepStatus;

  const G4VParticleChange* fParticleChange;
  const G4Track* fTrack;
  const G4TrackVector* fSecondary;
  const G4Step* fStep;
  G4StepPoint* fPreStepPoint;
  G4StepPoint* fPostStepPoint;

  const G4VPhysicalVolume* fCurrentVolume;
//  G4VSensitiveDetector* fSensitive;
  const G4VITProcess* fCurrentProcess;
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

  size_t fAtRestDoItProcTriggered;
  size_t fPostStepDoItProcTriggered;

  G4int fN2ndariesAtRestDoIt;
  G4int fN2ndariesAlongStepDoIt;
  G4int fN2ndariesPostStepDoIt;
  // These are the numbers of secondaries generated by the process
  // just executed.

//  G4Navigator *fNavigator;

  G4int fVerboseLevel;

  typedef std::vector<G4int> G4SelectedAtRestDoItVector;
  typedef std::vector<G4int> G4SelectedAlongStepDoItVector;
  typedef std::vector<G4int> G4SelectedPostStepDoItVector;
  G4SelectedAtRestDoItVector* fSelectedAtRestDoItVector;
  G4SelectedPostStepDoItVector* fSelectedPostStepDoItVector;

  G4double fPreviousStepSize;

  G4TouchableHandle fTouchableHandle;

//  G4SteppingControl StepControlFlag;

  G4double physIntLength;
  G4ForceCondition fCondition;
  G4GPILSelection fGPILSelection;
  // Above three variables are for the method
  // DefinePhysicalStepLength(). To pass these information to
  // the method Verbose, they are kept at here. Need a more
  // elegant mechanism.
};

#endif /* SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4VITSTEPPINGVERBOSE_HH_ */
