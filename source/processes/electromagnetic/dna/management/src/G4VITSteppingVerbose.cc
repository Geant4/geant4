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
 * G4VITSteppingVerbose.cc
 *
 *  Created on: Jun 22, 2015
 *      Author: mkaramit
 */

#include <G4VITSteppingVerbose.hh>
#include <G4UIcmdWithAnInteger.hh>

//______________________________________________________________________________

G4VITSteppingVerbose::G4VITSteppingVerbose()
{
  fpStepProcessor = nullptr;
  fpState = nullptr;
  fpProcessGeneralInfo = nullptr;

  PhysicalStep = -1;
  fStepStatus = fUndefined;

  fParticleChange = nullptr;
  fTrack = nullptr;
  fSecondary = nullptr;
  fStep = nullptr;
  fPreStepPoint = nullptr;
  fPostStepPoint = nullptr;

  fCurrentVolume = nullptr;
  //  fSensitive = fpStepProcessor->GetfSensitive();
  fCurrentProcess = nullptr;

  fAtRestDoItVector = nullptr;
  fAlongStepDoItVector = nullptr;
  fPostStepDoItVector = nullptr;

  fAtRestGetPhysIntVector = nullptr;
  fAlongStepGetPhysIntVector = nullptr;
  fPostStepGetPhysIntVector = nullptr;

  MAXofAtRestLoops = 0;
  MAXofAlongStepLoops = 0;
  MAXofPostStepLoops = 0;

  fAtRestDoItProcTriggered = 0;
  fPostStepDoItProcTriggered = 0;

  fN2ndariesAtRestDoIt = 0;
  fN2ndariesAlongStepDoIt = 0;
  fN2ndariesPostStepDoIt = 0;

  //  fNavigator = fpStepProcessor->GetfNavigator();

  fVerboseLevel = 0;
  fpVerboseUI = new G4UIcmdWithAnInteger("/chem/tracking/verbose", this);

  fSelectedAtRestDoItVector = nullptr;
  fSelectedPostStepDoItVector = nullptr;

  fPreviousStepSize = 0.;

  fTouchableHandle = nullptr;

  //  StepControlFlag = fpStepProcessor->GetStepControlFlag();

  physIntLength = 0;
  fCondition = InActivated;
  fGPILSelection = NotCandidateForSelection;

}

//______________________________________________________________________________

G4VITSteppingVerbose::~G4VITSteppingVerbose()
{
  delete fpVerboseUI;
}

//______________________________________________________________________________

void G4VITSteppingVerbose::CopyState()
{

  if(fpState != nullptr) *fpState = *(fpStepProcessor->GetProcessorState());
  else
  {
    fpState = new G4ITStepProcessorState(*fpStepProcessor->GetProcessorState());
  }

  fpProcessGeneralInfo = fpStepProcessor->GetCurrentProcessInfo();

  PhysicalStep = fpStepProcessor->GetPhysIntLength();
  fStepStatus = fpState->fStepStatus;

  fParticleChange = fpStepProcessor->GetParticleChange();
  fTrack = fpStepProcessor->GetTrack();
  fSecondary = fpStepProcessor->GetSecondaries();
  fStep = fpStepProcessor->GetStep();
  fPreStepPoint = fStep->GetPreStepPoint();
  fPostStepPoint = fStep->GetPostStepPoint();

  fCurrentVolume = fpStepProcessor->GetCurrentVolume();
//  fSensitive = fpStepProcessor->GetfSensitive();
  fCurrentProcess = fpStepProcessor->GetCurrentProcess();

  fAtRestDoItVector = fpProcessGeneralInfo->fpAtRestDoItVector;
  fAlongStepDoItVector = fpProcessGeneralInfo->fpAlongStepDoItVector;
  fPostStepDoItVector = fpProcessGeneralInfo->fpPostStepDoItVector;

  fAtRestGetPhysIntVector = fpProcessGeneralInfo->fpAtRestGetPhysIntVector;
  fAlongStepGetPhysIntVector =
      fpProcessGeneralInfo->fpAlongStepGetPhysIntVector;
  fPostStepGetPhysIntVector = fpProcessGeneralInfo->fpPostStepGetPhysIntVector;

  MAXofAtRestLoops = fpProcessGeneralInfo->MAXofAtRestLoops;
  MAXofAlongStepLoops = fpProcessGeneralInfo->MAXofAlongStepLoops;
  MAXofPostStepLoops = fpProcessGeneralInfo->MAXofPostStepLoops;

  fAtRestDoItProcTriggered = fpStepProcessor->GetAtRestDoItProcTriggered();
  fPostStepDoItProcTriggered = fpStepProcessor->GetPostStepDoItProcTriggered();

  fN2ndariesAtRestDoIt = fpStepProcessor->GetN2ndariesAtRestDoIt();
  fN2ndariesAlongStepDoIt = fpStepProcessor->GetN2ndariesAlongStepDoIt();
  fN2ndariesPostStepDoIt = fpStepProcessor->GetN2ndariesPostStepDoIt();

//  fNavigator = fpStepProcessor->GetfNavigator();

  fSelectedAtRestDoItVector = &(fpState->fSelectedAtRestDoItVector);
  fSelectedPostStepDoItVector = &(fpState->fSelectedPostStepDoItVector);

  fPreviousStepSize = fpState->fPreviousStepSize;

  fTouchableHandle = fpState->fTouchableHandle;

//  StepControlFlag = fpStepProcessor->GetStepControlFlag();

  physIntLength = fpStepProcessor->GetPhysIntLength();
  fCondition = fpStepProcessor->GetCondition();
  fGPILSelection = fpStepProcessor->GetGPILSelection();
}

//______________________________________________________________________________

void G4VITSteppingVerbose::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fpVerboseUI)
  {
    fVerboseLevel = fpVerboseUI->GetNewIntValue(newValue);
  }
}

//______________________________________________________________________________

G4String G4VITSteppingVerbose::GetCurrentValue(G4UIcommand* command)
{
  return command->ConvertToString(fVerboseLevel);
}

//______________________________________________________________________________

void G4VITSteppingVerbose::TrackingStarted(G4Track*
#ifdef G4VERBOSE
                                           track
#endif
)
{
#ifdef G4VERBOSE
  if(fVerboseLevel > 0)
  {
    TrackBanner(track, "G4ITTrackingManager::StartTracking : ");
  }
#endif

}

//______________________________________________________________________________

void G4VITSteppingVerbose::TrackingEnded(G4Track*
#ifdef G4VERBOSE
                                         track
#endif
)
{
#ifdef G4VERBOSE
  if(fVerboseLevel > 0)
  {
    TrackBanner(track, "G4ITTrackingManager::EndTracking : ");
  }
#endif
}

//______________________________________________________________________________

void G4VITSteppingVerbose::TrackBanner(G4Track* track, const G4String& message)
{
  G4cout << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  if(!message.empty())
  {
    G4cout << message;
  }
  G4cout << " * G4Track Information: "
         << "   Particle : " << track->GetDefinition()->GetParticleName()
         << ","
         << "   Track ID : " << track->GetTrackID()
         << ","
         << "   Parent ID : " << track->GetParentID()
         << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << G4endl;
}
