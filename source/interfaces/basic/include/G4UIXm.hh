// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIXm.hh,v 1.7 1999-11-10 15:01:15 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4UIXm_h
#define G4UIXm_h 

#if defined(G4UI_BUILD_XM_SESSION) || defined(G4UI_USE_XM)

#include "g4std/map"

#include <X11/Intrinsic.h>

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"

class G4UIsession;

// Class description :
//
//  G4UIXm : class to handle a Motif interactive session.
// G4UIXm is the Motif version of G4UIterminal.
//
//  A command box is at disposal for entering/recalling Geant4 commands.
//  A menubar could be customized through the AddMenu, AddButton methods.
//  Note that there are corresponding Geant4 commands to add a 
// menus in the menubar and add buttons in a menu.
//  Ex : 
//    /interactor/addMenu   test Test
//    /interactor/addButton test Init /run/initialize
//    /interactor/addButton test "Set gun" "/control/execute gun.g4m"
//    /interactor/addButton test "Run one event" "/run/beamOn 1"
//
//  Command completion, by typing "tab" key, is available on the 
// command line.
//
// Class description - end :

class G4UIXm : public G4VBasicShell, public G4VInteractiveSession {
public: // With description
  G4UIXm(int,char**);
  // (argv, argc) or (0, NULL) had to be given.
  G4UIsession* SessionStart();
  // To enter interactive X loop ; waiting/executing command,...
  void AddMenu(const char*,const char*);
  // To add a pulldown menu in the menu bar. 
  // First argument is the name of the menu.
  // Second argument is the label of the cascade button.
  // Ex : AddMenu("my_menu","My menu")
  void AddButton(const char*,const char*,const char*);
  // To add a push button in a pulldown menu.
  // First argument is the name of the menu.
  // Second argument is the label of the button.
  // Third argument is the Geant4 command executed when the button is fired.
  // Ex : AddButton("my_menu","Run","/run/beamOn 1"); 
public:
  ~G4UIXm();
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(G4String);
  G4int ReceiveG4cout(G4String);
  G4int ReceiveG4cerr(G4String);
  G4String GetCommand(Widget);
private:
  void SecondaryLoop(G4String);
  G4bool GetHelpChoice(G4int&);
  void ExitHelp();
private:
  Widget form,shell,command,menuBar,text;
  G4std::map<Widget,G4String,less<Widget> > commands;
  static void commandEnteredCallback(Widget,XtPointer,XtPointer);
  static void keyHandler(Widget,XtPointer,XEvent*,Boolean*);
  G4bool fHelp;
  G4int fHelpChoice;
  static void ButtonCallback(Widget,XtPointer,XtPointer);
};

#endif

#endif

