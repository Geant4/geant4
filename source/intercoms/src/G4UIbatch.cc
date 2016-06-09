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
// $Id: G4UIbatch.cc,v 1.14 2006/06/29 19:08:35 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

#include "G4UIbatch.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

G4bool G4UIbatch::commandFailed = false;

G4UIbatch::G4UIbatch(const char* fileName,G4UIsession* prevSession) 
 : previousSession(prevSession), macroFileName(fileName),
   openFailed(false)
{
  UImanager = G4UImanager::GetUIpointer();
  UImanager->SetSession(this);
  macroFile.open((char*)fileName);
  if(macroFile.fail())
  {
    G4cerr << "macro file <" << fileName << "> could not open."
         << G4endl;
    openFailed = true;
    commandFailed = true;
  }
  else
  {
    commandFailed = false;
  }
}

G4UIbatch::~G4UIbatch() 
{
  if(!openFailed) macroFile.close();
}

G4UIsession * G4UIbatch::SessionStart() 
{
  if(!openFailed)
  {
    char commandLine[256];
    G4int lineLength = 255;

    while(1)
    {
      macroFile.getline( commandLine, lineLength );
      if( macroFile.bad() )
      {
        G4cout << "Cannot read " << macroFileName << "." << G4endl;
        commandFailed = true;
        break;
      }
      if( macroFile.eof() ) break;
      commandLine[lineLength] = '\0';
      G4String commandString = commandLine;
      G4String nC= commandString.strip(G4String::both);
      if( commandLine[0] == '#')
      { if(G4UImanager::GetUIpointer()->GetVerboseLevel()==2)
        { G4cout << commandLine << G4endl; }
      }
      else if( nC.length() == 0 )
      { continue; }
      else if(nC == "exit")
      { break; }
      else
      { 
        G4int rc = UImanager->ApplyCommand(commandLine);
        if(rc)
        {
          switch(rc) 
          {
          case fCommandNotFound:
            G4cerr << "***** COMMAND NOT FOUND <"
                   << commandLine << "> *****" << G4endl;
            break;
          case fIllegalApplicationState:
            G4cerr << "***** Illegal application state <"
                   << commandLine << "> *****" << G4endl;
            break;
          default:
            G4int pn = rc%100;
            G4cerr << "***** Illegal parameter (" << pn << ") <"
                   << commandLine << "> *****" << G4endl;
          }
          G4cerr << "***** Command ignored *****" << G4endl;
          commandFailed = true;
        }
        if(commandFailed) break;
      }
    }
  }
  return previousSession;
}

void G4UIbatch::PauseSessionStart(G4String Prompt) 
{
  G4cout << "Pause session <" << Prompt << "> start." << G4endl;
  SessionStart();
  G4cout << "Pause session <" << Prompt << "> Terminate." << G4endl;
}


