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
// $Id: MyEventActionMessenger.cc,v 1.3 2001-07-11 09:56:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyEventActionMessenger.hh"

#include "MyEventAction.hh"
#include "G4UIparameter.hh"
#include "globals.hh"

MyEventActionMessenger::MyEventActionMessenger(MyEventAction * myEA)
:myEventAction(myEA)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand("/event/Draw",this);
  command->SetGuidance("Set drawFlag to Draw an event.");
  command->SetGuidance("  default value is 0 (false).");
  param = new G4UIparameter("flag (1:true/0:fasle)",'i',true);
  param->SetDefaultValue(0);
  command->SetParameter(param);
  AddUIcommand(command);
}

void MyEventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command->GetCommandName() == "Draw" )
  {
    G4int vl;
    const char* t = newValues;
    G4std::istrstream is((char*)t);
    is >> vl;
    myEventAction->SetDrawFlag(vl!=0);
  }
}

