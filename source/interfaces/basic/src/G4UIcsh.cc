// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcsh.cc,v 1.2 2000-06-14 03:18:59 asaim Exp $
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
G4String G4UIcsh::GetCommandLine()
//////////////////////////////////
{
  MakePrompt();
  G4cout << promptString << G4std::flush;

  G4String newCommand;
  newCommand.readLine(G4cin, FALSE);
  if (!G4cin.good()) {
    G4cin.clear(); 
    newCommand= "exit";
    return newCommand;
  }

  // multi-line
  while( newCommand[newCommand.length()-1] == '_' ) {
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

