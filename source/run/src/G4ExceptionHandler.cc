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
// G4ExceptionHandler implementation
//
// Author: M.Asai - August 2002
// --------------------------------------------------------------------

#include <stdlib.h>

#include "G4ExceptionHandler.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4String.hh"
#include "G4ios.hh"

#include "G4EventManager.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManagerKernel.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

// --------------------------------------------------------------------
G4ExceptionHandler::G4ExceptionHandler() {}

// --------------------------------------------------------------------
G4ExceptionHandler::~G4ExceptionHandler() {}

// --------------------------------------------------------------------
G4bool G4ExceptionHandler::operator==(const G4ExceptionHandler& right) const
{
  return (this == &right);
}

// --------------------------------------------------------------------
G4bool G4ExceptionHandler::operator!=(const G4ExceptionHandler& right) const
{
  return (this != &right);
}

// --------------------------------------------------------------------
G4bool G4ExceptionHandler::Notify(const char* originOfException,
                                  const char* exceptionCode,
                                  G4ExceptionSeverity severity,
                                  const char* description)
{
  static const G4String es_banner =
    "\n-------- EEEE ------- G4Exception-START -------- EEEE -------\n";
  static const G4String ee_banner =
    "\n-------- EEEE -------- G4Exception-END --------- EEEE -------\n";
  static const G4String ws_banner =
    "\n-------- WWWW ------- G4Exception-START -------- WWWW -------\n";
  static const G4String we_banner =
    "\n-------- WWWW -------- G4Exception-END --------- WWWW -------\n";
  std::ostringstream message;
  message << "*** G4Exception : " << exceptionCode << G4endl
          << "      issued by : " << originOfException << G4endl << description
          << G4endl;
  G4bool abortionForCoreDump = false;
  G4ApplicationState aps = G4StateManager::GetStateManager()->GetCurrentState();
  switch(severity)
  {
    case FatalException:
      G4cerr << es_banner << message.str()
             << "*** Fatal Exception *** core dump ***" << G4endl;
      DumpTrackInfo();
      G4cerr << ee_banner << G4endl;
      abortionForCoreDump = true;
      break;
    case FatalErrorInArgument:
      G4cerr << es_banner << message.str()
             << "*** Fatal Error In Argument *** core dump ***" << G4endl;
      DumpTrackInfo();
      G4cerr << ee_banner << G4endl;
      abortionForCoreDump = true;
      break;
    case RunMustBeAborted:
      if(aps == G4State_GeomClosed || aps == G4State_EventProc)
      {
        G4cerr << es_banner << message.str() << "*** Run Must Be Aborted ***"
               << G4endl;
        DumpTrackInfo();
        G4cerr << ee_banner << G4endl;
        G4RunManager::GetRunManager()->AbortRun(false);
      }
      abortionForCoreDump = false;
      break;
    case EventMustBeAborted:
      if(aps == G4State_EventProc)
      {
        G4cerr << es_banner << message.str() << "*** Event Must Be Aborted ***"
               << G4endl;
        DumpTrackInfo();
        G4cerr << ee_banner << G4endl;
        G4RunManager::GetRunManager()->AbortEvent();
      }
      abortionForCoreDump = false;
      break;
    default:
      G4cout << ws_banner << message.str()
             << "*** This is just a warning message. ***" << we_banner
             << G4endl;
      abortionForCoreDump = false;
      break;
  }
  return abortionForCoreDump;
}

// --------------------------------------------------------------------
void G4ExceptionHandler::DumpTrackInfo()
{
  const G4Track *theTrack = nullptr;
  const G4Step *theStep = nullptr;
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_EventProc)
  {
    G4SteppingManager *steppingMgr = G4RunManagerKernel::GetRunManagerKernel()
                                         ->GetTrackingManager()
                                         ->GetSteppingManager();
    theTrack = steppingMgr->GetfTrack();
    theStep = steppingMgr->GetfStep();
  }

  if (theTrack == nullptr)
  {
    G4cerr << " **** Track information is not available at this moment"
           << G4endl;
  }
  else
  {
    G4cerr << "G4Track (" << theTrack
           << ") - track ID = " << theTrack->GetTrackID()
           << ", parent ID = " << theTrack->GetParentID() << G4endl;
    G4cerr << " Particle type : "
           << theTrack->GetParticleDefinition()->GetParticleName();
    if(theTrack->GetCreatorProcess())
    {
      G4cerr << " - creator process : "
             << theTrack->GetCreatorProcess()->GetProcessName()
             << ", creator model : " << theTrack->GetCreatorModelName()
             << G4endl;
    }
    else
    {
      G4cerr << " - creator process : not available" << G4endl;
    }
    G4cerr << " Kinetic energy : "
           << G4BestUnit(theTrack->GetKineticEnergy(), "Energy")
           << " - Momentum direction : " << theTrack->GetMomentumDirection()
           << G4endl;
  }

  if (theStep == nullptr)
  {
    G4cerr << " **** Step information is not available at this moment"
           << G4endl;
  }
  else
  {
    G4cerr << " Step length : "
           << G4BestUnit(theStep->GetStepLength(), "Length")
           << " - total energy deposit : "
           << G4BestUnit(theStep->GetTotalEnergyDeposit(), "Energy") << G4endl;
    G4cerr << " Pre-step point : " << theStep->GetPreStepPoint()->GetPosition();
    G4cerr << " - Physical volume : ";
    if(theStep->GetPreStepPoint()->GetPhysicalVolume())
    {
      G4cerr << theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      if(theStep->GetPreStepPoint()->GetMaterial())
      {
        G4cerr << " (" << theStep->GetPreStepPoint()->GetMaterial()->GetName()
               << ")";
      }
      else
      {
        G4cerr << " (material not available)";
      }
    }
    else
    {
      G4cerr << "not available";
    }
    G4cerr << G4endl;
    if(theStep->GetPreStepPoint()->GetProcessDefinedStep())
    {
      G4cerr
        << " - defined by : "
        << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName()
        << " - step status : " << theStep->GetPreStepPoint()->GetStepStatus()
        << G4endl;
    }
    else
    {
      G4cerr << " - defined by : not available" << G4endl;
    }
    G4cerr << " Post-step point : "
           << theStep->GetPostStepPoint()->GetPosition();
    G4cerr << " - Physical volume : ";
    if(theStep->GetPostStepPoint()->GetPhysicalVolume())
    {
      G4cerr << theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      if(theStep->GetPostStepPoint()->GetMaterial())
      {
        G4cerr << " (" << theStep->GetPostStepPoint()->GetMaterial()->GetName()
               << ")";
      }
      else
      {
        G4cerr << " (material not available)";
      }
    }
    else
    {
      G4cerr << "not available";
    }
    G4cerr << G4endl;
    if(theStep->GetPostStepPoint()->GetProcessDefinedStep())
    {
      G4cerr << " - defined by : "
             << theStep->GetPostStepPoint()
                  ->GetProcessDefinedStep()
                  ->GetProcessName()
             << " - step status : "
             << theStep->GetPostStepPoint()->GetStepStatus() << G4endl;
    }
    else
    {
      G4cerr << " - defined by : not available" << G4endl;
    }
    G4cerr << " *** Note: Step information might not be properly updated."
           << G4endl;
  }
}
