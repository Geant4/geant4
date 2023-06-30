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

#if ! (defined(WIN32) || defined(__MINGW32__))

#  include "G4UIcommand.hh"
#  include "G4UIcommandTree.hh"
#  include "G4VUIshell.hh"

#  include <termios.h>

#  include <vector>

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

class G4UItcsh : public G4VUIshell
{
 public:
  G4UItcsh(const G4String& prompt = "%s> ", G4int maxhist = 100);
  ~G4UItcsh() override;

  void SetLsColor(TermColorIndex dirColor, TermColorIndex cmdColor) override;
  G4String GetCommandLineString(const char* msg = nullptr) override;

  void ResetTerminal() override;

 protected:
  void MakePrompt(const char* msg = nullptr) override;

  // Is cursor position at the last of command line ?
  G4bool IsCursorLast() const;

  void InitializeCommandLine();
  G4String ReadLine();

  // insert character
  void InsertCharacter(char cc);

  // backspace character
  void BackspaceCharacter();

  // delete character
  void DeleteCharacter();

  // clear command line
  void ClearLine();

  // clear after the cursor
  void ClearAfterCursor();

  // clear screen
  void ClearScreen();

  // move cursor forward
  void ForwardCursor();

  // move cursor backward
  void BackwardCursor();

  // move cursor to the top
  void MoveCursorTop();

  // move cursor to the end
  void MoveCursorEnd();

  // next command
  void NextCommand();

  // previous command
  void PreviousCommand();

  // list matched commands
  void ListMatchedCommand();

  // complete command
  void CompleteCommand();

  // utilities...
  G4String GetFirstMatchedString(const G4String& str1, const G4String& str2) const;

  void StoreHistory(G4String aCommand);

  // index is global history#
  G4String RestoreHistory(G4int index);

  void SetTermToInputMode();
  void RestoreTerm();

  G4String commandLine;  // command line string;
  G4int cursorPosition;  // cursor position
  G4String commandLineBuf;  // temp. command line;
  // history functionality  (history# is managed in itself)
  std::vector<G4String> commandHistory;
  G4int maxHistory;  // max# of histories stored
  G4int currentHistoryNo;  // global
  G4int relativeHistoryIndex;  // local index relative to current history#
  // (re)set termios
  termios tios;  // terminal mode (prestatus)
  G4String clearString;  // "clear code (^L)"
};

// ====================================================================
//   inline functions
// ====================================================================
inline G4bool G4UItcsh::IsCursorLast() const
{
  return cursorPosition == G4int(commandLine.length() + 1);
}

inline void G4UItcsh::SetLsColor(TermColorIndex dirColor, TermColorIndex cmdColor)
{
  lsColorFlag = true;
  directoryColor = dirColor;
  commandColor = cmdColor;
}

#endif
#endif
