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
// $Id: G4UIcsh.cc,v 1.6 2002-02-26 02:09:08 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4UIcsh.hh"

////////////////////////////////////////
G4UIcsh::G4UIcsh(const G4String& prompt)
 : G4VUIshell(prompt)
////////////////////////////////////////
{
}

///////////////////
G4UIcsh::~G4UIcsh()
///////////////////
{
}


//////////////////////////////////
G4String G4UIcsh::GetCommandLine(const char* msg)
//////////////////////////////////
{
  MakePrompt(msg);
  G4cout << promptString << G4std::flush;

  G4String newCommand;
  newCommand.readLine(G4cin, FALSE);
  if (!G4cin.good()) {
    G4cin.clear(); 
    newCommand= "exit";
    return newCommand;
  }
  newCommand = newCommand.strip(1,'\r'); // fix for odd behavior on Windows

  // multi-line
  while( (newCommand.length() > 0) &&
	 (newCommand[newCommand.length()-1] == '_') ) {
    G4String newLine;
    newCommand.remove(newCommand.length()-1);
    newLine.readLine(G4cin, FALSE);
    if (!G4cin.good()) { 
      G4cin.clear(); 
      newCommand= "exit";
      return newCommand;
    }
    newCommand.append(newLine);
  }
  
  return newCommand;
}

