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
#ifndef G4UIXm_h
#define G4UIXm_h 

#include <map>

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
//    /gui/addMenu   test Test
//    /gui/addButton test Init /run/initialize
//    /gui/addButton test "Set gun" "/control/execute gun.g4m"
//    /gui/addButton test "Run one event" "/run/beamOn 1"
//
//  Command completion, by typing "tab" key, is available on the 
// command line.
//
// Class description - end :

class G4UIXm : public G4VBasicShell, public G4VInteractiveSession {
public: // With description
  G4UIXm(G4int,char**);
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
  virtual void PauseSessionStart(const G4String&);
  virtual G4int ReceiveG4cout(const G4String&);
  virtual G4int ReceiveG4cerr(const G4String&);
  G4String GetCommand(Widget);
private:
  void SecondaryLoop(G4String);
  G4bool GetHelpChoice(G4int&);
  void ExitHelp() const;
  static void CommandEnteredCallback(Widget,XtPointer,XtPointer);
  static void keyHandler(Widget,XtPointer,XEvent*,Boolean*);
  static void ButtonCallback(Widget,XtPointer,XtPointer);

  Widget form,shell,command,menuBar,text;
  std::map<Widget,G4String, std::less<Widget> > commands;
  G4bool fHelp;
  G4int fHelpChoice;
  G4String menu_str[6] = { "form", "menuBar", "command",
                           "Clear", "clearButton", "text" };
};

#endif
