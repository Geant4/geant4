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
// $Id: G4VUIshell.hh 66892 2013-01-17 10:57:59Z gunter $
//

#ifndef G4VUIshell_h
#define G4VUIshell_h 1

#include "globals.hh"

// ====================================================================
//   Description: 
//   This class is the abstract base class for various UI shells.
//
//   GetCommadLineString() (virtual) returns a command string input from
//    a commad line.
//
//   Two pre-inplemented shell commands(still virtual) are also included, 
//   (somewhat differnt flavor from ones provided by G4VBasicShell)
//     ShowCurrentDirectory()   ... show current directory
//     ListCommand()            ... list commands
//
//   [prompt string substitution] (default)
//   %s ... current application status
//   %/ ... current working directory
//
// ====================================================================

// terminal color index
enum TermColorIndex{ BLACK=0, RED, GREEN, YELLOW, 
                     BLUE, PURPLE, CYAN, WHITE};

class G4UIcommandTree;

class G4VUIshell {
protected:
  G4String promptSetting; // including %-directive
  G4String promptString;
  virtual void MakePrompt(const char* msg=0);  // make prompt string
  G4int nColumn;  // column size of terminal (default=80)

  // color code support (effective if your terminal supports color code.)
  // default setting is off.
  G4bool lsColorFlag; // color flag for list command
  TermColorIndex directoryColor;
  TermColorIndex commandColor;

  // for treating G4 command tree...
  G4String currentCommandDir; // current command directory (absolute path)
  // get tree node
  G4UIcommandTree* GetCommandTree(const G4String& dir) const;  
  // absolute path name (ignore command)
  G4String GetAbsCommandDirPath(const G4String& apath) const;   
  // tail of path ( xxx/xxx/zzz -> zzz, trancated //// -> /)
  G4String GetCommandPathTail(const G4String& apath) const;  

public:
  G4VUIshell(const G4String& prompt="> ");
  virtual ~G4VUIshell();

  void SetNColumn(G4int ncol);
  void SetPrompt(const G4String& prompt);
  void SetCurrentDirectory(const G4String& ccd);
  virtual void SetLsColor(TermColorIndex, TermColorIndex);

  // shell commands
  virtual void ShowCurrentDirectory() const;
  virtual void ListCommand(const G4String& input, 
			   const G4String& candidate="") const;
  //  "candidate" is specified with full path.

  // get command string from a command line
  virtual G4String GetCommandLineString(const char* msg=0)= 0;

  virtual void ResetTerminal();
};

// ====================================================================
//   inline functions
// ====================================================================
inline void G4VUIshell::SetNColumn(G4int ncol)
{
  nColumn= ncol;
}

inline void G4VUIshell::SetPrompt(const G4String& prompt)
{
  promptSetting= prompt;
}

inline void G4VUIshell::SetCurrentDirectory(const G4String& dir)
{
  currentCommandDir= dir;
}

inline void G4VUIshell::SetLsColor(TermColorIndex, TermColorIndex)
{
}

inline void G4VUIshell::ShowCurrentDirectory() const
{
  G4cout << currentCommandDir << G4endl;
}

#endif
