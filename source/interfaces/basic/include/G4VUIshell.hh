// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUIshell.hh,v 1.2 2000-07-22 10:52:28 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VUIshell_h
#define G4VUIshell_h 1

#include "globals.hh"

//
//   Description: 
//   This class is the abstract base class for various UI shells.
//
//   GetCommadLine() (virtual) returns a command string input from
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

// terminal color index
enum TermColorIndex{ BLACK=0, RED, GREEN, YELLOW, 
                     BLUE, PURPLE, CYAN, WHITE};

class G4UIcommandTree;

class G4VUIshell {
protected:
  G4String promptString;
  G4String promptSetting; // including %-directive
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
  ~G4VUIshell();

  void SetNColumn(G4int ncol);
  void SetPrompt(const G4String& prompt);
  void SetCurrentDirectory(const G4String& ccd);

  // shell commands
  virtual void ShowCurrentDirectory() const;
  virtual void ListCommand(const G4String& input, 
			   const G4String& candidate="") const;
  //  "candidate" is specified with full path.

  // get command string from a command line
  virtual G4String GetCommandLine(const char* msg=0)= 0;
};


// inlines...
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

inline void G4VUIshell::ShowCurrentDirectory() const
{
  G4cout << currentCommandDir << G4endl;
}

#endif
