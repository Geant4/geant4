// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIWin32.hh,v 1.7 1999-11-10 15:20:44 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4UIWin32_h
#define G4UIWin32_h 

#if defined(G4UI_BUILD_WIN32_SESSION) || defined(G4UI_USE_WIN32)

#include "g4std/map"
#include "g4std/vector"

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"

#include "G4Win32.hh"

class G4VInteractorManager;

// Class description :
//
//  G4UIWin32 : class to handle a Windows interactive session.
// G4UIWin32 is the Windows version of G4UIterminal.
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

class G4UIWin32 : public G4VBasicShell, public G4VInteractiveSession {
public: // With description
  G4UIWin32 (HINSTANCE,HINSTANCE,LPSTR,int);
  // The WinMain arguments have to be given.
  G4UIsession* SessionStart      ();
  // To enter interactive Win32 loop ; waiting/executing command,...
  void AddMenu (const char*,const char*);
  // To add a pulldown menu in the menu bar. 
  // First argument is the name of the menu.
  // Second argument is the label of the cascade button.
  // Ex : AddMenu("my_menu","My menu")
  void AddButton (const char*,const char*,const char*);
  // To add a push button in a pulldown menu.
  // First argument is the name of the menu.
  // Second argument is the label of the button.
  // Third argument is the Geant4 command executed when the button is fired.
  // Ex : AddButton("my_menu","Run","/run/beamOn 1"); 
public:
  ~G4UIWin32 ();
  void Prompt (G4String);
  void SessionTerminate ();
  void PauseSessionStart (G4String);
  G4int ReceiveG4cout(G4String);
  G4int ReceiveG4cerr(G4String);
  G4String GetCommand (int);
  void TextAppendString(char*);
private:
  void SecondaryLoop (G4String);
  G4bool GetHelpChoice(G4int&);
  void ExitHelp();
private:
  G4VInteractorManager* interactorManager;
  HWND mainWindow;
  HWND textWindow,editWindow;
  HMENU menuBar,defaultMenu;
  G4std::map<int,G4String, G4std::less<int> > commands;
  void* textBuffer;
  int textRows,textCols;
  static LRESULT CALLBACK MainWindowProc(HWND,UINT,WPARAM,LPARAM);
  static LRESULT CALLBACK TextWindowProc(HWND,UINT,WPARAM,LPARAM);
  static LRESULT CALLBACK EditWindowProc(HWND,UINT,WPARAM,LPARAM);
  G4bool fHelp;
  G4int fHelpChoice;
  G4std::vector<G4String> fHistory;
  int fHistoryPos;
};

#endif

#endif

