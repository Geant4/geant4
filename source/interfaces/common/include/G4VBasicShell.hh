// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VBasicShell.hh,v 1.5 1999-12-15 14:50:47 gunter Exp $
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

