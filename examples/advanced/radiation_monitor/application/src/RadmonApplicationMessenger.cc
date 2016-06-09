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
// File name:     RadmonApplicationMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationMessenger.cc,v 1.1.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-04 $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/application/"

// Include files
#include "RadmonApplicationMessenger.hh"
#include "RadmonApplicationEventNumbering.hh"
#include "RadmonApplicationEventTracks.hh"
#include "RadmonApplicationRunNumbering.hh"
#include "RadmonEventAction.hh"
#include "RadmonRunAction.hh"
#include "Randomize.hh"



                                                RadmonApplicationMessenger :: RadmonApplicationMessenger()
:
 RadmonMessenger(COMMANDS_PATH, "Interactive application options."),
 eventNumbering(new  RadmonApplicationEventNumbering),
 eventTracks(new RadmonApplicationEventTracks),
 runNumbering(new RadmonApplicationRunNumbering),
 RADMON_INITIALIZE_COMMAND(EnableRunsDump),
 RADMON_INITIALIZE_COMMAND(DisableRunsDump),
 RADMON_INITIALIZE_COMMAND(DumpEventsEvery),
 RADMON_INITIALIZE_COMMAND(DisableEventsDump),
 RADMON_INITIALIZE_COMMAND(EnableTracksVisualisation),
 RADMON_INITIALIZE_COMMAND(DisableTracksVisualisation),
 RADMON_INITIALIZE_COMMAND(SetSeed)
{
 RadmonEventAction * eventAction(RadmonEventAction::Instance());
 eventAction->AttachObserver(eventNumbering);
 eventAction->AttachObserver(eventTracks);
 
 RadmonRunAction * runAction(RadmonRunAction::Instance());
 runAction->AttachObserver(runNumbering);

 RADMON_CREATE_COMMAND_0ARGS(EnableRunsDump,              "Disables the run number dump");
 RADMON_CREATE_COMMAND_0ARGS(DisableRunsDump,             "Enables the run number dump");
 RADMON_CREATE_COMMAND_1ARG (DumpEventsEvery,             "Enables the event number dump every n events", "nEvents");
 RADMON_CREATE_COMMAND_0ARGS(DisableEventsDump,           "Disables the event number dump");

 #ifdef    G4VIS_USE
  RADMON_CREATE_COMMAND_0ARGS(EnableTracksVisualisation,  "Enables the tracks plotting at the end of the event");
  RADMON_CREATE_COMMAND_0ARGS(DisableTracksVisualisation, "Disables the tracks plotting at the end of the event");
 #endif /* G4VIS_USE */

 RADMON_CREATE_COMMAND_1ARG (SetSeed,                     "Sets the seed of the random generator", "seed");
}



                                                RadmonApplicationMessenger :: ~RadmonApplicationMessenger()
{
 RADMON_DESTROY_COMMAND(SetSeed);

 #ifdef    G4VIS_USE
  RADMON_DESTROY_COMMAND(DisableTracksVisualisation);
  RADMON_DESTROY_COMMAND(EnableTracksVisualisation);
 #endif /* G4VIS_USE */

 RADMON_DESTROY_COMMAND(DisableEventsDump);
 RADMON_DESTROY_COMMAND(DumpEventsEvery);
 RADMON_DESTROY_COMMAND(DisableRunsDump);
 RADMON_DESTROY_COMMAND(EnableRunsDump);

 RadmonEventAction * eventAction(RadmonEventAction::Instance());
 eventAction->DetachObserver(eventNumbering);
 delete eventNumbering;
 eventAction->DetachObserver(eventTracks);
 delete eventTracks;
 
 RadmonRunAction * runAction(RadmonRunAction::Instance());
 runAction->DetachObserver(runNumbering);
 delete runNumbering;
}





G4String                                        RadmonApplicationMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonApplicationMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonApplicationMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(EnableRunsDump)
  RADMON_SET_COMMAND(DisableRunsDump)
  RADMON_SET_COMMAND(DumpEventsEvery)
  RADMON_SET_COMMAND(DisableEventsDump)
  RADMON_SET_COMMAND(EnableTracksVisualisation)
  RADMON_SET_COMMAND(DisableTracksVisualisation)
  RADMON_SET_COMMAND(SetSeed)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonApplicationMessenger :: OnEnableRunsDump(const G4String & /* value */)
{
 runNumbering->Enable();
}



void                                            RadmonApplicationMessenger :: OnDisableRunsDump(const G4String & /* value */)
{
 runNumbering->Disable();
}



void                                            RadmonApplicationMessenger :: OnDumpEventsEvery(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, & args))
  return; 

 G4int every(G4UIcommand::ConvertToInt(args));
 
 if (every<=0)
 {
  G4cout << "RadmonApplicationMessenger::OnDumpEventsEvery: nEvents must be positive." << G4endl;
  return;
 }

 eventNumbering->SetDumpEvery(every);
}



void                                            RadmonApplicationMessenger :: OnDisableEventsDump(const G4String & /* value */)
{
 eventNumbering->SetDumpEvery(0);
}



void                                            RadmonApplicationMessenger :: OnEnableTracksVisualisation(const G4String & /* value */)
{
 eventTracks->Enable();
}



void                                            RadmonApplicationMessenger :: OnDisableTracksVisualisation(const G4String & /* value */)
{
 eventTracks->Disable();
}



void                                            RadmonApplicationMessenger :: OnSetSeed(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, & args))
  return; 

 CLHEP::HepRandom::setTheSeed(G4UIcommand::ConvertToInt(args));
}
