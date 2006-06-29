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
// $Id: G4VBasicShell.hh,v 1.7 2006-06-29 19:10:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VBasicShell_H
#define G4VBasicShell_H 1

class G4UIcommandTree;
class G4UIcommand;

#include "G4UIsession.hh"
#include "globals.hh"

// Class description :
//
//  G4VBasicShell : a base class to extract common things to various
// sessions. 
//
//  It handles "seek" completion logic, help logic.
//  VBasicShell handles also commands like "cd, ls, pwd"
// without passing by a Geant4 "intercom" command. This feature,
// which is similar to a UNIX shell one, had given.
// its name to the class.
//
// Class description - end :

class G4VBasicShell : public G4UIsession
{
  public:
    G4VBasicShell();
    virtual ~G4VBasicShell();

    virtual G4UIsession* SessionStart() = 0;
    // null should be returned for interactive session

    virtual void PauseSessionStart(G4String Prompt) = 0;
    // Prompt string can be ignored

  protected:
    G4String ModifyToFullPathCommand(const char* aCommandLine);
    // convert "BeamOn 10" to "/run/BeamOn 10" if the
    // current working directory is "/run/"

    G4String GetCurrentWorkingDirectory();
    // directory string starts with '/' and ends with '/'

    G4bool ChangeDirectory(const char* newDir);
    // change directory to newDir
    // false will be returned if the target directory doesn't exist

    G4UIcommandTree* FindDirectory(const char* dirName);
    // find G4UIcommandTree object
    // null returned if the taregt does not exist

    G4UIcommand* FindCommand(const char* commandName);
    // find G4UIcommand object
    // null returned if the target does not exist

    G4String Complete(G4String);
    // command completion
    G4String FindMatchingPath(G4UIcommandTree*,G4String);

    /////////////////////////////////////////////
    // Methods involving an interactive G4cout //
    /////////////////////////////////////////////
    virtual void ExecuteCommand(G4String);
    virtual G4bool GetHelpChoice(G4int&) = 0;
    virtual void ExitHelp() = 0;
    void ApplyShellCommand(G4String,G4bool&,G4bool&);
    void ShowCurrent(G4String);
    void ChangeDirectoryCommand(G4String);
    void ListDirectory(G4String);
    void TerminalHelp(G4String);

  private:
    G4String currentDirectory;

    G4String ModifyPath(G4String tempPath);
};

#endif

