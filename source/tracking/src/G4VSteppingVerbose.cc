// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSteppingVerbose.cc,v 1.5 2001-02-08 07:39:54 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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

G4VSteppingVerbose* G4VSteppingVerbose::fInstance = 0;
G4VSteppingVerbose::G4VSteppingVerbose(){;}
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

   currentMinimumStep = fManager->GetcurrentMinimumStep();
   numberOfInteractionLengthLeft = fManager->GetnumberOfInteractionLengthLeft();

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

   fTouchable1 = fManager->GetfTouchable1();
   fTouchable2 = fManager->GetfTouchable2();
   fIsTouchable1Free = fManager->GetfIsTouchable1Free();
   fIsTouchable2Free = fManager->GetfIsTouchable2Free();

   StepControlFlag = fManager->GetStepControlFlag();

   physIntLength = fManager->GetphysIntLength();
   fCondition = fManager->GetfCondition();
   fGPILSelection = fManager->GetfGPILSelection();
}


