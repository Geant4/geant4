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


///////////////////////////////////////////////////////
G4String G4UIcsh::GetCommandLineString(const char* msg)
///////////////////////////////////////////////////////
{
  MakePrompt(msg);
  G4cout << promptString << std::flush;

  G4String newCommand;
  G4StrUtil::readline(G4cin, newCommand, false);
  if (!G4cin.good()) {
    G4cin.clear(); 
    newCommand= "exit";
    return newCommand;
  }
  G4StrUtil::rstrip(newCommand, '\r'); // fix for odd behavior on Windows

  // multi-line
  while( (newCommand.length() > 0) &&
         (newCommand[G4int(newCommand.length()-1)] == '_') ) {
    G4String newLine;
    newCommand.erase(newCommand.length()-1);
    G4StrUtil::readline(G4cin, newLine, false);
    if (!G4cin.good()) { 
      G4cin.clear(); 
      newCommand= "exit";
      return newCommand;
    }
    newCommand.append(newLine);
  }
  
  return newCommand;
}
