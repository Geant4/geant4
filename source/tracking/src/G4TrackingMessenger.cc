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
// $Id: G4TrackingMessenger.cc 77241 2013-11-22 09:55:47Z gcosmo $
//
//---------------------------------------------------------------
//
// G4TrackingMessenger.cc
//
// Description:
//   This is a messenger class to interface to exchange information
//   between tracking and stepping managers and UI.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "G4TrackingMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4TrackStatus.hh"
#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4IdentityTrajectoryFilter.hh"

///////////////////////////////////////////////////////////////////
G4TrackingMessenger::G4TrackingMessenger(G4TrackingManager * trMan)
///////////////////////////////////////////////////////////////////
: trackingManager(trMan)
{
  steppingManager = trackingManager->GetSteppingManager();

  TrackingDirectory = new G4UIdirectory("/tracking/");
  TrackingDirectory->SetGuidance("TrackingManager and SteppingManager control commands.");

  AbortCmd = new G4UIcmdWithoutParameter("/tracking/abort",this);
  AbortCmd->SetGuidance("Abort current G4Track processing.");


  ResumeCmd = new G4UIcmdWithoutParameter("/tracking/resume",this);
  ResumeCmd->SetGuidance("Resume current G4Track processing.");

  StoreTrajectoryCmd = new G4UIcmdWithAnInteger("/tracking/storeTrajectory",this);
  StoreTrajectoryCmd->SetGuidance("Store trajectories or not.");
  StoreTrajectoryCmd->SetGuidance(" 0 : Don't Store trajectories.");
  StoreTrajectoryCmd->SetGuidance(" !=0 : Store trajectories.");
  StoreTrajectoryCmd->SetGuidance(" 1 : Choose G4Trajectory as default.");
  StoreTrajectoryCmd->SetGuidance(" 2 : Choose G4SmoothTrajectory as default.");
  StoreTrajectoryCmd->SetGuidance(" 3 : Choose G4RichTrajectory as default.");
  StoreTrajectoryCmd->SetGuidance(" 4 : Choose G4RichTrajectory with auxiliary points as default.");
  StoreTrajectoryCmd->SetParameterName("Store",true);
  StoreTrajectoryCmd->SetDefaultValue(0);
  StoreTrajectoryCmd->SetRange("Store >=0 && Store <= 4"); 


  VerboseCmd = new G4UIcmdWithAnInteger("/tracking/verbose",this);
#ifdef G4VERBOSE
  VerboseCmd->SetGuidance("Set Verbose level of tracking category.");
  VerboseCmd->SetGuidance(" -1 : Silent.");
  VerboseCmd->SetGuidance(" 0 : Silent.");
  VerboseCmd->SetGuidance(" 1 : Minium information of each Step.");
  VerboseCmd->SetGuidance(" 2 : Addition to Level=1, info of secondary particles.");
  VerboseCmd->SetGuidance(" 3 : Addition to Level=1, pre/postStepoint information");
  VerboseCmd->SetGuidance("     after all AlongStep/PostStep process executions.");
  VerboseCmd->SetGuidance(" 4 : Addition to Level=3, pre/postStepoint information");
  VerboseCmd->SetGuidance("     at each AlongStepPostStep process execuation."); 
  VerboseCmd->SetGuidance(" 5 : Addition to Level=4, proposed Step length information");
  VerboseCmd->SetGuidance("     from each AlongStepPostStep process."); 
  VerboseCmd->SetParameterName("verbose_level",true);
  VerboseCmd->SetDefaultValue(0);
  VerboseCmd->SetRange("verbose_level >=-1  ");
#else 
  VerboseCmd->SetGuidance("You need to recompile the tracking category defining G4VERBOSE ");  
#endif
}

////////////////////////////////////////////
G4TrackingMessenger::~G4TrackingMessenger()
////////////////////////////////////////////
{
  delete TrackingDirectory;
  delete AbortCmd;
  delete ResumeCmd;
  delete StoreTrajectoryCmd;
  delete VerboseCmd;
}

///////////////////////////////////////////////////////////////////////////////
void G4TrackingMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
///////////////////////////////////////////////////////////////////////////////
{
  if( command == VerboseCmd ){
    trackingManager->SetVerboseLevel(VerboseCmd->ConvertToInt(newValues));
  }

  if( command  == AbortCmd ){
    steppingManager->GetTrack()->SetTrackStatus(fStopAndKill);
    G4UImanager::GetUIpointer()->ApplyCommand("/control/exit");
  }

  if( command  == ResumeCmd ){
        G4UImanager::GetUIpointer()->ApplyCommand("/control/exit");
  }

  static G4ThreadLocal G4IdentityTrajectoryFilter* auxiliaryPointsFilter = 0;
  if(!auxiliaryPointsFilter) auxiliaryPointsFilter = new G4IdentityTrajectoryFilter;
  if( command == StoreTrajectoryCmd ){
    G4int trajType = StoreTrajectoryCmd->ConvertToInt(newValues);
    if(trajType==2||trajType==4)
    {
      G4TransportationManager::GetTransportationManager()
        ->GetPropagatorInField()->SetTrajectoryFilter(auxiliaryPointsFilter);
    }
    else
    {
      G4TransportationManager::GetTransportationManager()
        ->GetPropagatorInField()->SetTrajectoryFilter(0);
    }
    trackingManager->SetStoreTrajectory(trajType);
  }
}


////////////////////////////////////////////////////////////////////
G4String G4TrackingMessenger::GetCurrentValue(G4UIcommand * command)
////////////////////////////////////////////////////////////////////
{
  if( command == VerboseCmd ){
    return VerboseCmd->ConvertToString(trackingManager->GetVerboseLevel());
  }
  else if( command == StoreTrajectoryCmd ){
    return StoreTrajectoryCmd->ConvertToString(trackingManager->GetStoreTrajectory());
  }
  return G4String('\0');
}













