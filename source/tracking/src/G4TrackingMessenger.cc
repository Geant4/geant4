// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackingMessenger.cc,v 1.1 1999-01-07 16:14:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

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
  StoreTrajectoryCmd->SetGuidance(" 1 : Store trajectories.");
  StoreTrajectoryCmd->SetGuidance(" 0 : Don't Store trajectories.");
  StoreTrajectoryCmd->SetParameterName("Store",true);
  StoreTrajectoryCmd->SetDefaultValue(0);
  StoreTrajectoryCmd->SetRange("Store >=0 && Store <= 1"); 


  VerboseCmd = new G4UIcmdWithAnInteger("/tracking/verbose",this);
#ifdef G4VERBOSE
  VerboseCmd->SetGuidance("Set Verbose level of tracking category.");
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
  VerboseCmd->SetRange("verbose_level >=0  ");
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
    G4int vl;
    const char* t = newValues;
    istrstream is((char*)t);
    is >> vl;
    trackingManager->SetVerboseLevel(vl);
  }

  if( command  == AbortCmd ){
    steppingManager->GetTrack()->SetTrackStatus(fStopAndKill);
    G4UImanager::GetUIpointer()->ApplyCommand("/control/exit");
  }

  if( command  == ResumeCmd ){
        G4UImanager::GetUIpointer()->ApplyCommand("/control/exit");
  }

  if( command == StoreTrajectoryCmd ){
    G4int vl;
    const char* t = newValues;
    istrstream is((char*)t);
    is >> vl;
    trackingManager->SetStoreTrajectory(vl!=0);
  }
}


////////////////////////////////////////////////////////////////////
G4String G4TrackingMessenger::GetCurrentValue(G4UIcommand * command)
////////////////////////////////////////////////////////////////////
{
  if( command == VerboseCmd ){
    char line[100];
    ostrstream os(line,100);
    os << trackingManager->GetVerboseLevel() << '\0';
    return G4String(line);
  }
  return G4String('\0');
}













