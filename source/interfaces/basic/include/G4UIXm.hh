// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIXm.hh,v 1.2 1999-04-13 01:26:25 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4UIXm_h
#define G4UIXm_h 

#if defined(G4UI_BUILD_XM_SESSION) || defined(G4UI_USE_XM)

#include <rw/tvhdict.h>

#include <X11/Intrinsic.h>

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"

class G4UIsession;

class G4UIXm : public G4VBasicShell, public G4VInteractiveSession {
public:
  G4UIXm(int,char**);
  ~G4UIXm();
  G4UIsession* SessionStart();
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(G4String);
  G4int ReceiveG4cout(G4String);
  G4int ReceiveG4cerr(G4String);
  void AddMenu(const char*,const char*);
  void AddButton(const char*,const char*,const char*);
  G4String GetCommand(Widget);
private:
  void SecondaryLoop(G4String);
  G4bool GetHelpChoice(G4int&);
  void ExitHelp();
private:
  Widget shell,command,menuBar,text;
  RWTValHashDictionary<Widget,G4String> commands;
  static void commandCallback(Widget,XtPointer,XtPointer);
  G4bool fHelp;
  G4int fHelpChoice;
  static void ButtonCallback(Widget,XtPointer,XtPointer);
};

#endif

#endif

