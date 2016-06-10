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
// $Id: G4VSteppingVerbose.cc 66872 2013-01-15 01:25:57Z japost $
//
//---------------------------------------------------------------
//
// G4VSteppingVerbose.cc
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "G4VSteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"

G4ThreadLocal G4VSteppingVerbose* G4VSteppingVerbose::fInstance = 0;
G4ThreadLocal G4int G4VSteppingVerbose::Silent = 0;
G4ThreadLocal G4int G4VSteppingVerbose::SilentStepInfo = 0;

G4VSteppingVerbose::G4VSteppingVerbose()
  : fManager(0), fUserSteppingAction(0), 
    PhysicalStep(0.), GeometricalStep(0.), CorrectedStep(0.),
    PreStepPointIsGeom(false), FirstStep(false),
    TempInitVelocity(0.), TempVelocity(0.), Mass(0.), sumEnergyChange(0.),
    fParticleChange(0), fTrack(0), fSecondary(0), fStep(0), 
    fPreStepPoint(0), fPostStepPoint(0), fCurrentVolume(0),
    fSensitive(0), fCurrentProcess(0), fAtRestDoItVector(0),
    fAlongStepDoItVector(0), fPostStepDoItVector(0), fAtRestGetPhysIntVector(0),
    fAlongStepGetPhysIntVector(0), fPostStepGetPhysIntVector(0),
    MAXofAtRestLoops(0), MAXofAlongStepLoops(0), MAXofPostStepLoops(0),
    currentMinimumStep(0.), numberOfInteractionLengthLeft(0.),
    fAtRestDoItProcTriggered(0), fAlongStepDoItProcTriggered(0),
    fPostStepDoItProcTriggered(0), fN2ndariesAtRestDoIt(0),
    fN2ndariesAlongStepDoIt(0), fN2ndariesPostStepDoIt(0),
    fNavigator(0), verboseLevel(0), fSelectedAtRestDoItVector(0),
    fSelectedAlongStepDoItVector(0), fSelectedPostStepDoItVector(0),
    fPreviousStepSize(0.), physIntLength(0.)
{
  if(fInstance!= 0)
  {
    G4Exception("G4VSteppingVerbose::G4VSteppingVerbose()",
                "Tracking0014", FatalException,
                "Only one SteppingVerbose class can be instantiated.");
  }
}
G4VSteppingVerbose::~G4VSteppingVerbose(){;}

//////////////////////////////////////////////////////////////////
void G4VSteppingVerbose::SetManager(G4SteppingManager* const fMan)
//////////////////////////////////////////////////////////////////
{
  fManager=fMan;
}

//////////////////////////////////////////////////
void G4VSteppingVerbose::CopyState()
//////////////////////////////////////////////////
{

   fUserSteppingAction = fManager->GetUserAction();
   //   fVerbose = this;

   PhysicalStep = fManager->GetPhysicalStep();
   GeometricalStep = fManager->GetGeometricalStep();
   CorrectedStep = fManager->GetCorrectedStep();
   PreStepPointIsGeom = fManager->GetPreStepPointIsGeom();
   FirstStep = fManager->GetFirstStep();
   fStepStatus = fManager->GetfStepStatus();

   TempInitVelocity = fManager->GetTempInitVelocity();
   TempVelocity = fManager->GetTempVelocity();
   Mass = fManager->GetMass();

   sumEnergyChange = fManager->GetsumEnergyChange();

   fParticleChange = fManager->GetfParticleChange();
   fTrack = fManager->GetfTrack(); 
   fSecondary = fManager->GetfSecondary();
   fStep = fManager->GetfStep();
   fPreStepPoint = fManager->GetfPreStepPoint();
   fPostStepPoint = fManager->GetfPostStepPoint();

   fCurrentVolume = fManager->GetfCurrentVolume();
   fSensitive = fManager->GetfSensitive();
   fCurrentProcess = fManager->GetfCurrentProcess();

   fAtRestDoItVector = fManager->GetfAtRestDoItVector(); 
   fAlongStepDoItVector = fManager->GetfAlongStepDoItVector();
   fPostStepDoItVector = fManager->GetfPostStepDoItVector();

   fAtRestGetPhysIntVector = fManager->GetfAtRestGetPhysIntVector();
   fAlongStepGetPhysIntVector = fManager->GetfAlongStepGetPhysIntVector();
   fPostStepGetPhysIntVector = fManager->GetfPostStepGetPhysIntVector();

   MAXofAtRestLoops = fManager->GetMAXofAtRestLoops();
   MAXofAlongStepLoops = fManager->GetMAXofAlongStepLoops();
   MAXofPostStepLoops = fManager->GetMAXofPostStepLoops();

   fAtRestDoItProcTriggered = fManager->GetfAtRestDoItProcTriggered();
   fAlongStepDoItProcTriggered = fManager->GetfAlongStepDoItProcTriggered();
   fPostStepDoItProcTriggered = fManager->GetfPostStepDoItProcTriggered();

   fN2ndariesAtRestDoIt = fManager->GetfN2ndariesAtRestDoIt();
   fN2ndariesAlongStepDoIt = fManager->GetfN2ndariesAlongStepDoIt();
   fN2ndariesPostStepDoIt = fManager->GetfN2ndariesPostStepDoIt();

   fNavigator = fManager->GetfNavigator();

   verboseLevel = fManager->GetverboseLevel();

   fSelectedAtRestDoItVector = fManager->GetfSelectedAtRestDoItVector();
   fSelectedAlongStepDoItVector = fManager->GetfSelectedAlongStepDoItVector();
   fSelectedPostStepDoItVector = fManager->GetfSelectedPostStepDoItVector();

   fPreviousStepSize = fManager->GetfPreviousStepSize();

   fTouchableHandle = fManager->GetTouchableHandle();

   StepControlFlag = fManager->GetStepControlFlag();

   physIntLength = fManager->GetphysIntLength();
   fCondition = fManager->GetfCondition();
   fGPILSelection = fManager->GetfGPILSelection();
}

void G4VSteppingVerbose::SetInstance(G4VSteppingVerbose* Instance)
{
  fInstance=Instance;
}

G4VSteppingVerbose* G4VSteppingVerbose::GetInstance()
{
  return fInstance;
}

G4int G4VSteppingVerbose::GetSilent()
{
  return Silent;
}

void G4VSteppingVerbose::SetSilent(G4int fSilent)
{
  Silent=fSilent;
}

G4int G4VSteppingVerbose::GetSilentStepInfo()
{
  return SilentStepInfo;
}

void G4VSteppingVerbose::SetSilentStepInfo(G4int fSilent)
{
  SilentStepInfo=fSilent;
}

