// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIxvtMessenger.cc,v 1.2 1999-04-13 01:25:03 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifdef G4UI_BUILD_XVT_SESSION

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The original code (written by Makoto Asai) has been adapted for use      //
// with the XVT Geant GUI                                                   //
//                                                                          //
// By: Simon Prior 11/08/97                                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "G4UIxvtMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIxvt.hh"
#include "G4ios.hh"


G4UIxvtMessenger::G4UIxvtMessenger(G4UIxvt* tempxvt)
:xvtptr(tempxvt)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand("/control/setBreakAt",this);
  command->SetGuidance("Set a break point at the begining or end of");
  command->SetGuidance("each run or event.");
  command->SetGuidance(" Available points :");
  command->SetGuidance("  beginOfRun, endOfRun, beginOfEvent, EndOfEvent");
  param = new G4UIparameter("point",'c',false);
  command->SetParameter(param);
  AddUIcommand( command );

  command = new G4UIcommand("/control/removeBreakPoint",this);
  command->SetGuidance("Reset all of assigned break.");
  AddUIcommand( command );

  command = new G4UIcommand("/control/noticeState",this);
  command->SetGuidance("Toggle Verbose flag of G4UIxvt for state change.");
  command->SetGuidance("If it is on, all state change will be displayed.");
  AddUIcommand( command );
}

void G4UIxvtMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if(command->GetCommandName()=="setBreakAt")
  {
    int id = -1;
    if(newValue=="beginOfRun") id=0;
    if(newValue=="beginOfEvent") id=1;
    if(newValue=="EndOfEvent") id=2;
    if(newValue=="endOfRun") id=3;
    if(id<0)
    { G4cout << "Unknown break point <" << newValue << "> ignored." << endl; }
    else
    { xvtptr->set_breakPointAt(id,true); }
  }
  if(command->GetCommandName()=="removeBreakPoint")
  {
    for(int i=0;i<4;i++)
    { xvtptr->set_breakPointAt(i,false); }
  }
  if(command->GetCommandName()=="noticeState")
  {
    xvtptr->set_verbose();
  }
}

#endif
