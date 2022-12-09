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

#ifndef G4UItcsh_h
#define G4UItcsh_h 1

#if !(defined(WIN32) || defined(__MINGW32__))

#include <termios.h>
#include <vector>
#include "G4VUIshell.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"

// ====================================================================
//   Description:
//   This class gives tcsh-like shell.
//   
//   If your terminal supports color code, colored strings are available 
//   in ListCommand(). For activating color support, 
//   e.g.
//     tcsh-> SetLsColor(GREEN, CYAN); // (dir, command) color
//   
//   [key binding]
//   ^A ... move cursor to the top
//   ^B ... backward cursor ([LEFT])
//   ^D ... delete/exit/show matched list
//   ^E ... move cursor to the end
//   ^F ... forward cursor ([RIGHT])
//   ^K ... clear after the cursor
//   ^L ... clear screen (not implemented)
//   ^N ... next command ([DOWN])
//   ^P ... previous command ([UP])
//   TAB... command completion
//   DEL... backspace
//   BS ... backspace
//
//   [prompt string substitution]
//   %s ... current application status
//   %/ ... current working directory
//   %h ... history# (different from G4 history#)
//
// ====================================================================

class G4UItcsh : public G4VUIshell {
protected:
  virtual void MakePrompt(const char* msg=0);

  G4String commandLine;    // command line string;
  G4int cursorPosition;    // cursor position 
  G4String commandLineBuf; // temp. command line;
  G4bool IsCursorLast() const; 
                           // Is cursor position at the last of command line ?

  void InitializeCommandLine();
  G4String ReadLine();
  void InsertCharacter(char cc); // insert character
  void BackspaceCharacter();     // backspace character
  void DeleteCharacter();        // delete character
  void ClearLine();              // clear command line
  void ClearAfterCursor();       // clear after the cursor
  void ClearScreen();            // clear screen

  void ForwardCursor();          // move cursor forward
  void BackwardCursor();         // move cursor backward
  void MoveCursorTop();          // move cursor to the top
  void MoveCursorEnd();          // move cursor to the end

  void NextCommand();            // next command
  void PreviousCommand();        // previous command

  void ListMatchedCommand();     // list matched commands
  void CompleteCommand();        // complete command  
  
  // utilities...
  G4String GetFirstMatchedString(const G4String& str1, 
				 const G4String& str2) const;

  // history functionality  (history# is managed in itself)
  std::vector<G4String> commandHistory;
  G4int maxHistory;            // max# of histories stored
  G4int currentHistoryNo;      // global
  G4int relativeHistoryIndex;  // local index relative to current history#

  void StoreHistory(G4String aCommand);
  G4String RestoreHistory(G4int index); // index is global history#


  // (re)set termios
  termios tios; // terminal mode (prestatus)
  G4String clearString;  // "clear code (^L)"
  void SetTermToInputMode();
  void RestoreTerm();

public:
  G4UItcsh(const G4String& prompt="%s> ", G4int maxhist=100);
  ~G4UItcsh();
  
  void SetLsColor(TermColorIndex dirColor, TermColorIndex cmdColor);
  virtual G4String GetCommandLineString(const char* msg=0);

  virtual void ResetTerminal();
};

// ====================================================================
//   inline functions
// ====================================================================
inline G4bool G4UItcsh::IsCursorLast() const
{
  if(cursorPosition == G4int(commandLine.length()+1)) return TRUE;
  else return FALSE;
}

inline void G4UItcsh::SetLsColor(TermColorIndex dirColor, 
				 TermColorIndex cmdColor)
{
  lsColorFlag= TRUE;
  directoryColor= dirColor;
  commandColor= cmdColor;
}

#endif
#endif

