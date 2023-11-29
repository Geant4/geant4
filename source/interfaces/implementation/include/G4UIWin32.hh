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
// G4UIWin32
//
// Class description :
//
// G4UIWin32 is a class to handle a Windows interactive session.
// It is the Windows version of G4UIterminal.
//
//  A command box is at disposal for entering/recalling Geant4 commands.
//  A menubar could be customized through the AddMenu, AddButton methods.
//  Note that there are corresponding Geant4 commands to add a
//  menu in the menubar and add buttons in a menu.
//  Ex :
//    /gui/addMenu   test Test
//    /gui/addButton test Init /run/initialize
//    /gui/addButton test "Set gun" "/control/execute gun.g4m"
//    /gui/addButton test "Run one event" "/run/beamOn 1"
//
// Command completion, by typing "tab" key, is available on the
// command line.

// Original author: G.Barrand, 1998
// Revised: O.Pena-Rodriguez, March 2021
// -------------------------------------------------------------------
#ifndef G4UIWin32_h
#define G4UIWin32_h

#define ID_EXIT_APP 40001
#define ID_OPEN_MACRO 40002
#define ID_SAVE_VIEWER_STATE 40003
#define ID_RUN_BEAMON 40004
#define ID_RUN_CMD 40005
#define ID_VIEW_SOLID 40006
#define ID_VIEW_WIREFRAME 40007
#define ID_PROJ_ORTHOGRAPHIC 40008
#define ID_PROJ_PERSPECTIVE 40009
#define ID_ZOOM_IN 40010
#define ID_ZOOM_OUT 40011
#define ID_ORIENTATION_XY 40012
#define ID_ORIENTATION_XZ 40013
#define ID_ORIENTATION_YZ 40014
#define ID_ORIENTATION_OBLIQUE 40015
#define ID_HELP_ABOUT 40016
#define ID_LOG_CLEAN 40017
#define ID_LOG_SAVE 40018

#define NUM_BUTTONS 12
#define MAX_HISTORY_ITEMS 10

#define IDC_MAIN_EDIT 101
#define IDC_MAIN_TOOL 102
#define IDC_MAIN_TREE_VIEW 104
#define IDC_MAIN_COMBO 105
#define IDC_MAIN_STATUS 106

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"
#include "G4Win32.hh"

#include <CommCtrl.h>

#include <map>
#include <vector>

class G4VInteractorManager;

class G4UIWin32 : public G4VBasicShell, public G4VInteractiveSession
{
 public:
  G4UIWin32();
  ~G4UIWin32() override;

  // To enter interactive Win32 loop ; waiting/executing command,...
  G4UIsession* SessionStart() override;

  // To add a pulldown menu in the menu bar.
  // First argument is the name of the menu.
  // Second argument is the label of the cascade button.
  // Ex : AddMenu("my_menu","My menu")
  void AddMenu(const char*, const char*) override;

  // To add a push button in a pulldown menu.
  // First argument is the name of the menu.
  // Second argument is the label of the button.
  // Third argument is the Geant4 command executed when the button is fired.
  // Ex : AddButton("my_menu","Run","/run/beamOn 1");
  void AddButton(const char*, const char*, const char*) override;

  void Prompt(const G4String&);
  void SessionTerminate();
  void PauseSessionStart(const G4String&) override;

  G4int ReceiveG4debug(const G4String&) override;
  G4int ReceiveG4cout(const G4String&) override;
  G4int ReceiveG4cerr(const G4String&) override;

  G4String GetCommand(G4int);
  // void TextAppendString(char*);

 private:
  void SecondaryLoop(const G4String&);
  G4bool GetHelpChoice(G4int&) override;
  void ExitHelp() const override;

  G4bool CreateComponents(HWND);
  G4bool ResizeComponents(HWND);
  void ProcessTabKey();
  void ProcessEscKey();
  void ProcessEnterKey();
  void ProcessUpKey();
  void ProcessDownKey();

  G4bool ProcessDefaultCommands(G4int);
  static G4String GetToolTips(G4int);
  G4String GetHelpTreeToolTips(HTREEITEM);

  static G4String ConvertNewLines(G4String);

  void HelpTreeDoubleClick(HTREEITEM);

  G4bool SaveLogFile(LPCTSTR);
  void AddText(LPSTR);

  void DoOpenMacro(HWND);
  void DoSaveViewer(HWND);
  void DoSaveLog(HWND);

  G4bool InitHelpTreeItems();
  HTREEITEM AddItemToHelpTree(LPTSTR, HTREEITEM = TVI_ROOT);
  static G4String GetShortCommandPath(G4String);
  LPSTR GetItemPath(HTREEITEM);

  void CreateHelpTree(HTREEITEM, G4UIcommandTree*);

 private:
  HWND fHWndMainWindow;
  HWND fHWndEditor;
  HWND fHWndToolBar;
  HWND fHWndComboBox;
  HWND fHWndComboEditor;
  HWND fHWndHelpTree;
  HWND fHWndStatus;

  G4VInteractorManager* interactorManager;
  HMENU menuBar;
  std::map<G4int, G4String, std::less<G4int>> commands;

  static LRESULT CALLBACK MainWindowProc(HWND, UINT, WPARAM, LPARAM);
  // New wndproc for the combo box
  static LRESULT CALLBACK ComboEditorWindowProc(HWND, UINT, WPARAM, LPARAM);

  G4bool fHelp;
  G4int fHelpChoice;
  std::vector<G4String> fHistory;
  G4int fHistoryPos;
};

#endif
