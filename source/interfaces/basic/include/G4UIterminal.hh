// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIterminal.hh,v 1.6 2000-07-22 10:52:28 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIterminal_h
#define G4UIterminal_h 1

#include "g4std/fstream"
#include "G4UImanager.hh"
#include "G4VBasicShell.hh"
#include "G4VUIshell.hh"

//
//   Description:
//
//   This class inherits the class G4UIsession.
//   This is the class to use a character-terminal sesion.
//
//   Usage:  
//       G4UIsession* terminalSession = new G4UIterminal; 
//   or  G4UIsession* terminalSession = new G4UIterminal(new your-shell); 
//
//     A character-terminal session  "terminalSession" is instantiated.
//     G4cout stream is redirected by default to the constructed instance.
//
//   terminalSession-> SessionStart(); // "terminalSession" is started.
//   delete terminalSession;           // "terminalSession"  is deleted.
//
//
//   In default(no arguments are given), csh-like shell is instantiated.
//   If you want to use another shell (e.g. tcsh-like), you can give
//   your favorite shell in an argument of the constructor.
//
//   Which shell? / How to define your own shell?
//   Currently two kinds of shells,
//                   G4UIcsh / G4UItcsh
//   , are presented. 
//   They inherit the abstract base class, G4VUIshell.
//   In order to define your own shell, 
//     - Define your own shell class derived from G4VUIshell.
//     - Implement GetCommandLine() method (pure virtual).
//     - Add more functionality, if need.
//
//   For more detail, see source codes.
//

class G4UIterminal : public G4VBasicShell {
private:
  G4UImanager* UI;
  // shell
  G4VUIshell* shell;

  // program states
  G4bool iExit;
  G4bool iCont;

public:
  G4UIterminal(G4VUIshell* aShell=0);
  ~G4UIterminal();

  void SetPrompt(const G4String& prompt);

private:
  void ExecuteCommand(G4String aCommand);
  G4String GetCommand(const char* msg=0);
  G4bool GetHelpChoice(G4int& aInt);
  void ExitHelp();
  
public:
  // These methods are implementation of corresponding virtual methods
  // of G4UIsession class.
  virtual G4UIsession* SessionStart();  
  virtual void PauseSessionStart(G4String msg);
  virtual G4int ReceiveG4cout(G4String coutString);
  virtual G4int ReceiveG4cerr(G4String cerrString);

};

#endif

