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
// $Id: G4UIterminal.hh 66892 2013-01-17 10:57:59Z gunter $
//
// ====================================================================
//   G4UIterminal.cc
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
// ====================================================================
#ifndef G4UIterminal_h
#define G4UIterminal_h 1

#include <fstream>
#include "G4UImanager.hh"
#include "G4VBasicShell.hh"
#include "G4VUIshell.hh"

class G4UIterminal : public G4VBasicShell {
private:
  G4UImanager* UI;
  // shell
  G4VUIshell* shell;

  // program states
  G4bool iExit;
  G4bool iCont;

public:
  G4UIterminal(G4VUIshell* aShell=0, G4bool qsig=true);
  ~G4UIterminal();

  void SetPrompt(const G4String& prompt);

private:
  virtual void ExecuteCommand(const G4String& aCommand);
  virtual G4bool GetHelpChoice(G4int& aInt);
  virtual void ExitHelp() const;
  G4String GetCommand(const char* msg=0);
  
public:
  // These methods are implementation of corresponding virtual methods
  // of G4UIsession class.
  virtual G4UIsession* SessionStart();  
  virtual void PauseSessionStart(const G4String& msg);
  virtual G4int ReceiveG4cout(const G4String& coutString);
  virtual G4int ReceiveG4cerr(const G4String& cerrString);
};

#endif

