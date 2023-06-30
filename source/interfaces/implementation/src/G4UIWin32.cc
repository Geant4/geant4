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
// Original author: G.Barrand, 1998
// Rewrite: O.Pena-Rodriguez, March 2021
// --------------------------------------------------------------------

#include "G4UIWin32.hh"

#include "G4MTcoutDestination.hh"
#include "G4StateManager.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UImanager.hh"
#include "G4Win32.hh"

#include <commctrl.h>
#include <windows.h>
#include <wingdi.h>

#include <cstring>
#include <utility>

/***************************************************************************/
static char mainClassName[] = "G4UIWin32";
static G4bool exitSession = true;
static G4bool exitPause = true;
static G4bool exitHelp = true;
static G4UIsession* tmpSession = nullptr;

// Original wndproc for the Combo Editor
static WNDPROC origComboEditorWindowProc;

static G4bool ConvertStringToInt(const char*, G4int&);

static G4int actionIdentifier = 0;

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4UIWin32::G4UIWin32()
  : fHWndMainWindow(nullptr),
    fHWndEditor(nullptr),
    fHWndToolBar(nullptr),
    fHWndComboBox(nullptr),
    fHWndComboEditor(nullptr),
    fHWndHelpTree(nullptr),
    fHWndStatus(nullptr),
    menuBar(nullptr),
    fHelp(false),
    fHelpChoice(0),
    fHistoryPos(-1)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI != nullptr) {
    UI->SetSession(this);
    UI->SetG4UIWindow(this);
  }

  // Ensure that the common control DLL is loaded
  INITCOMMONCONTROLSEX commCtrls;
  commCtrls.dwSize = sizeof(INITCOMMONCONTROLSEX);
  commCtrls.dwICC = ICC_BAR_CLASSES | ICC_LISTVIEW_CLASSES;
  // this loads list-view and toolbar classes.
  InitCommonControlsEx(&commCtrls);

  interactorManager = G4Win32::getInstance();

  static G4bool Done = false;
  if (! Done) {
    WNDCLASS wc;
    wc.style = 0;
    wc.lpfnWndProc = (WNDPROC)G4UIWin32::MainWindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = ::GetModuleHandle(nullptr);
    wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.lpszMenuName = mainClassName;
    wc.lpszClassName = mainClassName;

    if (! RegisterClass(&wc)) {
      MessageBox(nullptr, "G4UIWin32: Win32 window registration failed!", "Error!",
        MB_ICONEXCLAMATION | MB_OK);
      G4cout << "G4UIWin32: Win32 window registration failed!" << G4endl;
      return;
    }

    Done = true;
  }

  menuBar = CreateMenu();

  // Add some initial options to the menu
  HMENU hMenu = CreatePopupMenu();
  AppendMenu(menuBar, MF_POPUP, (UINT_PTR)hMenu, "&Geant4");
  AddInteractor("Geant4", (G4Interactor)hMenu);

  AppendMenu(hMenu, MF_STRING, ID_OPEN_MACRO, "&Open macro...");
  AppendMenu(hMenu, MF_STRING, ID_SAVE_VIEWER_STATE, "&Save viewer state...");
  AppendMenu(hMenu, MF_SEPARATOR, -1, "");
  AppendMenu(hMenu, MF_STRING, ID_RUN_BEAMON, "&Beam On");
  AppendMenu(hMenu, MF_SEPARATOR, -1, "");
  AppendMenu(hMenu, MF_STRING, ID_EXIT_APP, "E&xit");

  hMenu = CreatePopupMenu();
  AppendMenu(menuBar, MF_POPUP, (UINT_PTR)hMenu, "&View");
  AddInteractor("View", (G4Interactor)hMenu);

  AppendMenu(hMenu, MF_STRING, ID_VIEW_SOLID, "S&olid");
  AppendMenu(hMenu, MF_STRING, ID_VIEW_WIREFRAME, "&Wireframe");
  AppendMenu(hMenu, MF_SEPARATOR, -1, "");
  AppendMenu(hMenu, MF_STRING, ID_PROJ_ORTHOGRAPHIC, "&Orthographic");
  AppendMenu(hMenu, MF_STRING, ID_PROJ_PERSPECTIVE, "P&erspective");
  AppendMenu(hMenu, MF_SEPARATOR, -1, "");
  AppendMenu(hMenu, MF_STRING, ID_ORIENTATION_XY, "&X-Y Plane");
  AppendMenu(hMenu, MF_STRING, ID_ORIENTATION_XZ, "X-&Z Plane");
  AppendMenu(hMenu, MF_STRING, ID_ORIENTATION_YZ, "&Y-Z Plane");
  AppendMenu(hMenu, MF_STRING, ID_ORIENTATION_OBLIQUE, "&Oblique");

  hMenu = CreatePopupMenu();
  AppendMenu(menuBar, MF_POPUP, (UINT_PTR)hMenu, "&Zoom");
  AddInteractor("Zoom", (G4Interactor)hMenu);

  AppendMenu(hMenu, MF_STRING, ID_ZOOM_IN, "Zoom &In");
  AppendMenu(hMenu, MF_STRING, ID_ZOOM_OUT, "Zoom &Out");

  tmpSession = this;
  fHWndMainWindow = ::CreateWindowEx(WS_EX_CLIENTEDGE, mainClassName, "Geant4",
    WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
    CW_USEDEFAULT, nullptr, menuBar, ::GetModuleHandle(nullptr), nullptr);

  if (fHWndMainWindow == nullptr) {
    MessageBox(nullptr, "Window Creation Failed!", "Error!", MB_ICONEXCLAMATION | MB_OK);
    return;
  }
  tmpSession = nullptr;
  ::SetWindowLongPtr(fHWndMainWindow, GWLP_USERDATA, (LONG_PTR)this);

  ::SetForegroundWindow(fHWndMainWindow);
  ::ShowWindow(fHWndMainWindow, SW_SHOWDEFAULT);
  ::UpdateWindow(fHWndMainWindow);

  if (UI != nullptr) UI->SetCoutDestination(this);

  // TODO: Manage multithreaded output
  // #ifdef G4MULTITHREADED
  // explicitly request that cout/cerr messages from threads are ALSO propagated to the master.
  // masterG4coutDestination = this;
  // #endif
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4UIWin32::~G4UIWin32()
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI != nullptr) {
    UI->SetSession(nullptr);
    UI->SetG4UIWindow(nullptr);
    UI->SetCoutDestination(nullptr);
  }

  // TODO: Manage multithreaded output
  // #ifdef G4MULTITHREADED
  // masterG4coutDestination = 0; // set to cout when UI is deleted
  // #endif

  if (fHWndStatus != nullptr) ::SetWindowLongPtr(fHWndStatus, GWLP_USERDATA, LONG(NULL));
  if (fHWndHelpTree != nullptr) ::SetWindowLongPtr(fHWndHelpTree, GWLP_USERDATA, LONG(NULL));
  if (fHWndComboBox != nullptr) ::SetWindowLongPtr(fHWndComboBox, GWLP_USERDATA, LONG(NULL));
  if (fHWndToolBar != nullptr) ::SetWindowLongPtr(fHWndToolBar, GWLP_USERDATA, LONG(NULL));
  if (fHWndEditor != nullptr) ::SetWindowLongPtr(fHWndEditor, GWLP_USERDATA, LONG(NULL));
  if (fHWndMainWindow != nullptr) {
    ::SetWindowLongPtr(fHWndMainWindow, GWLP_USERDATA, LONG(NULL));
    ::DestroyWindow(fHWndMainWindow);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4UIsession* G4UIWin32::SessionStart()
{
  if (interactorManager != nullptr) {
    Prompt("session");
    exitSession = false;

    // TODO: Ensure that the list of commands is updated
    // Load commands into Help Tree View
    InitHelpTreeItems();

    interactorManager->DisableSecondaryLoop();
    void* event;
    while ((event = interactorManager->GetEvent()) != nullptr) {
      interactorManager->DispatchEvent(event);
      if (exitSession) break;
    }
    interactorManager->EnableSecondaryLoop();
    return this;
  }
  else
    return this;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::Prompt(const G4String& a_prompt) {}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::SessionTerminate() {}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::PauseSessionStart(const G4String& a_state)
{
  if (a_state == "G4_pause> ") {
    SecondaryLoop("Pause, type continue to exit this state");
  }

  if (a_state == "EndOfEvent") {
    // Picking with feed back in event data Done here !!!
    SecondaryLoop("End of event, type continue to exit this state");
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::SecondaryLoop(const G4String& a_prompt)
{
  if (interactorManager != nullptr) {
    Prompt(a_prompt);
    exitPause = false;
    void* event;
    while ((event = interactorManager->GetEvent()) != nullptr) {
      interactorManager->DispatchEvent(event);
      if (exitPause) break;
    }
    Prompt("session");
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4int G4UIWin32::ReceiveG4debug(const G4String& a_string)
{
  // Geant4 uses UNIX's style for new lines (\n)
  // we must convert them to Windows' style (\r\n)
  G4String str = ConvertNewLines(a_string);

  AddText((LPSTR)str.data());

  return 0;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4int G4UIWin32::ReceiveG4cout(const G4String& a_string)
{
  // Geant4 uses UNIX's style for new lines (\n)
  // we must convert them to Windows' style (\r\n)
  G4String str = ConvertNewLines(a_string);

  AddText((LPSTR)str.data());

  return 0;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4int G4UIWin32::ReceiveG4cerr(const G4String& a_string)
{
  // Geant4 uses UNIX's style for new lines (\n)
  // we must convert them to Windows' style (\r\n)
  G4String str = ConvertNewLines(a_string);

  AddText((LPSTR)str.data());

  return 0;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool G4UIWin32::GetHelpChoice(G4int& aInt)
{
  fHelp = true;

  if (interactorManager != nullptr) {
    Prompt("Help");
    exitHelp = false;
    void* event;
    while ((event = interactorManager->GetEvent()) != nullptr) {
      interactorManager->DispatchEvent(event);
      if (exitHelp) break;
    }
    Prompt("session");
    //
    if (! fHelp) return false;
    aInt = fHelpChoice;
    fHelp = false;
    return true;
  }
  else
    return false;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::ExitHelp() const {}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::AddMenu(const char* a_name, const char* a_label)
{
  if (a_name != nullptr) {
    HMENU hMenu = CreatePopupMenu();
    AppendMenu(menuBar, MF_POPUP, (UINT_PTR)hMenu, a_label);
    AddInteractor(a_name, (G4Interactor)hMenu);
    DrawMenuBar(fHWndMainWindow);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::AddButton(const char* a_menu, const char* a_label, const char* a_command)
{
  if ((a_menu != nullptr) && (a_label != nullptr) && (a_command != nullptr)) {
    HMENU hMenu = (HMENU)GetInteractor(a_menu);
    actionIdentifier++;
    commands[actionIdentifier] = a_command;
    AppendMenu(hMenu, MF_STRING, actionIdentifier, a_label);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4String G4UIWin32::GetCommand(G4int a_id) { return commands[a_id]; }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
LRESULT CALLBACK G4UIWin32::MainWindowProc(
  HWND aWindow, UINT aMessage, WPARAM wParam, LPARAM lParam)
{
  switch (aMessage) {
    case WM_CREATE: {
      auto* This = (G4UIWin32*)tmpSession;
      if (This != nullptr) {
        if (! This->CreateComponents(aWindow)) {
          MessageBox(aWindow, "Could not create components.", "Error", MB_OK | MB_ICONERROR);
          return false;
        }
      }
    }
      return 0;

    case WM_SIZE: {
      auto* This = (G4UIWin32*)::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      if (This != nullptr) {
        if (! This->ResizeComponents(aWindow)) {
          MessageBox(aWindow, "Could not resize components.", "Error", MB_OK | MB_ICONERROR);
          return false;
        }
      }
    }
      return 0;

    case WM_CLOSE:
      DestroyWindow(aWindow);
      return 0;

    case WM_DESTROY:
      PostQuitMessage(0);
      return 0;

    case WM_SETFOCUS: {
      auto* This = (G4UIWin32*)::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      if (This != nullptr) SetFocus(This->fHWndComboBox);
    }
      return 0;

    case WM_NOTIFY: {
      auto* This = (G4UIWin32*)::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      if (This != nullptr) {
        switch (((LPNMHDR)lParam)->code) {
          // Tooltip for Toolbar
          case TTN_NEEDTEXT: {
            auto lpttt = (LPTOOLTIPTEXT)lParam;
            lpttt->hinst = nullptr;
            UINT idButton = lpttt->hdr.idFrom;
            lpttt->lpszText = (LPSTR)This->GetToolTips(idButton).c_str();
          } break;

          // Tooltip for TreeView
          case TVN_GETINFOTIP: {
            auto pTip = (LPNMTVGETINFOTIP)lParam;
            pTip->pszText = (LPSTR)This->GetHelpTreeToolTips(pTip->hItem).c_str();
          } break;

          // Double click for TreeView
          case NM_DBLCLK: {
            auto lpnmh = (LPNMHDR)lParam;
            auto item = TreeView_GetSelection(lpnmh->hwndFrom);
            This->HelpTreeDoubleClick(item);
          } break;
        }
      }
    }
      return 0;

    case WM_COMMAND: {
      auto* This = (G4UIWin32*)::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      if (This != nullptr)
        if (! This->ProcessDefaultCommands(LOWORD(wParam)))
          // If the command was not processed, do it now
          switch (LOWORD(wParam)) {
            case IDC_MAIN_EDIT: {
              // We have to release some space when the buffer is full
              if (HIWORD(wParam) == EN_ERRSPACE || HIWORD(wParam) == EN_MAXTEXT) {
                G4int bufferSize =
                  SendMessage(This->fHWndEditor, EM_GETLIMITTEXT, (WPARAM)0, (LPARAM)0);

                // Select the first third of the text
                SendMessage(This->fHWndEditor, EM_SETSEL, (WPARAM)0, (LPARAM)bufferSize / 3);
                // Delete it
                SendMessage(This->fHWndEditor, EM_REPLACESEL, (WPARAM)0, (LPARAM) "");
                // Scroll to the bottom
                SendMessage(This->fHWndEditor, WM_VSCROLL, SB_BOTTOM, NULL);
              }
            } break;
            default:
              if (! This->fHelp) {
                G4String command = This->GetCommand(wParam);
                This->ApplyShellCommand(command, exitSession, exitPause);
              }
          }
    }
      return 0;
    default:
      //  For all the other cases, call the default window procedure.
      return DefWindowProc(aWindow, aMessage, wParam, lParam);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
LRESULT CALLBACK G4UIWin32::ComboEditorWindowProc(
  HWND aWindow, UINT aMessage, WPARAM wParam, LPARAM lParam)
{
  // We need to go two steps up: Editor -> ComboBox -> Window
  HWND parent = GetParent(GetParent(aWindow));
  auto* This = (G4UIWin32*)::GetWindowLongPtr(parent, GWLP_USERDATA);

  switch (aMessage) {
    case WM_KEYDOWN:
      switch (wParam) {
        case VK_TAB: {
          if (This != nullptr) {
            if (This->fHelp) break;

            This->ProcessTabKey();
          }
        }
          return 0;  // Do not jump into origComboEditorWindowProc.

        case VK_ESCAPE: {
          if (This != nullptr) This->ProcessEscKey();
        }
          return 0;  // Do not jump into origComboEditorWindowProc.

        case VK_RETURN: {
          if (This != nullptr) This->ProcessEnterKey();
        }
          return 0;  // Do not jump into origComboEditorWindowProc.

        case VK_UP: {
          if (This != nullptr) This->ProcessUpKey();
        }
          return 0;  // Do not jump into origComboEditorWindowProc.

        case VK_DOWN: {
          if (This != nullptr) This->ProcessDownKey();
        }
          return 0;  // Do not jump into origComboEditorWindowProc.
      }
      break;

    case WM_KEYUP:
    case WM_CHAR:
      switch (wParam) {
        case VK_TAB:
        case VK_ESCAPE:
        case VK_RETURN:
        case VK_UP:
        case VK_DOWN:
          return 0;  // Do not jump into origComboEditorWindowProc.
      }
  }

  //  Call the original window procedure for default processing.
  return CallWindowProc(origComboEditorWindowProc, aWindow, aMessage, wParam, lParam);
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool G4UIWin32::CreateComponents(HWND aWindow)
{
  HFONT hfDefault;
  TBBUTTON tbb[NUM_BUTTONS];
  TBADDBITMAP tbab;
  RECT rcClient;  // dimensions of client area

  G4int statwidths[] = {100, -1};

  // Create Edit Control
  fHWndEditor = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "",
    WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_HSCROLL | ES_MULTILINE | ES_AUTOVSCROLL |
      ES_AUTOHSCROLL | ES_READONLY,
    0, 0, 100, 100, aWindow, (HMENU)IDC_MAIN_EDIT, GetModuleHandle(nullptr), nullptr);
  if (fHWndEditor == nullptr)
    MessageBox(aWindow, "Could not create edit box.", "Error", MB_OK | MB_ICONERROR);

  // Set editor font
  // hfDefault = (HFONT) GetStockObject(DEFAULT_GUI_FONT);
  hfDefault = CreateFont(-10, -8, 0, 0, 0, false, 0, 0, OEM_CHARSET, OUT_RASTER_PRECIS,
    CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, FIXED_PITCH, "System");
  SendMessage(fHWndEditor, WM_SETFONT, (WPARAM)hfDefault, MAKELPARAM(false, 0));

  // Set editor's buffer size (the default value is too small)
  SendMessage(fHWndEditor, EM_SETLIMITTEXT, (WPARAM)500000, (LPARAM)0);

  // Create Toolbar
  fHWndToolBar = CreateWindowEx(0, TOOLBARCLASSNAME, nullptr,
    WS_CHILD | WS_VISIBLE | TBSTYLE_FLAT | TBSTYLE_TOOLTIPS, 0, 0, 0, 0, aWindow,
    (HMENU)IDC_MAIN_TOOL, GetModuleHandle(nullptr), nullptr);
  if (fHWndToolBar == nullptr)
    MessageBox(aWindow, "Could not create tool bar.", "Error", MB_OK | MB_ICONERROR);

  // Required for backward compatibility.
  SendMessage(fHWndToolBar, TB_BUTTONSTRUCTSIZE, (WPARAM)sizeof(TBBUTTON), (LPARAM)0);

  // Load standard images
  tbab.hInst = HINST_COMMCTRL;
  tbab.nID = IDB_STD_SMALL_COLOR;
  SendMessage(fHWndToolBar, TB_ADDBITMAP, (WPARAM)0, (LPARAM)&tbab);

  // Load history images
  tbab.hInst = HINST_COMMCTRL;
  tbab.nID = IDB_HIST_SMALL_COLOR;
  SendMessage(fHWndToolBar, TB_ADDBITMAP, (WPARAM)0, (LPARAM)&tbab);

  G4int btnBMP[NUM_BUTTONS] = {STD_FILEOPEN, STD_FILESAVE, -1, STD_FIND, STD_FIND, -1,
    15 + HIST_FORWARD, -1, STD_HELP, -1, STD_FILENEW, STD_FILESAVE};
  G4int btnSTL[NUM_BUTTONS] = {TBSTYLE_BUTTON, TBSTYLE_BUTTON, TBSTYLE_SEP, TBSTYLE_BUTTON,
    TBSTYLE_BUTTON, TBSTYLE_SEP, TBSTYLE_BUTTON, TBSTYLE_SEP, TBSTYLE_BUTTON, TBSTYLE_SEP,
    TBSTYLE_BUTTON, TBSTYLE_BUTTON};
  G4int btnCMD[NUM_BUTTONS] = {ID_OPEN_MACRO, ID_SAVE_VIEWER_STATE, -1, ID_ZOOM_IN, ID_ZOOM_OUT, -1,
    ID_RUN_BEAMON, -1, ID_HELP_ABOUT, -1, ID_LOG_CLEAN, ID_LOG_SAVE};
  ZeroMemory(tbb, sizeof(tbb));
  for (G4int i = 0; i < NUM_BUTTONS; i++) {
    tbb[i].iBitmap = btnBMP[i];
    tbb[i].fsState = TBSTATE_ENABLED;
    tbb[i].fsStyle = btnSTL[i];
    tbb[i].idCommand = btnCMD[i];
  }

  SendMessage(fHWndToolBar, TB_ADDBUTTONS, sizeof(tbb) / sizeof(TBBUTTON), (LPARAM)&tbb);

  // Create the Combobox
  fHWndComboBox = CreateWindowEx(0, WC_COMBOBOX, TEXT(""),
    CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, 150, 0, 200, 200,
    aWindow, (HMENU)IDC_MAIN_COMBO, GetModuleHandle(nullptr), nullptr);

  // Display an initial item in the selection field
  SendMessage(fHWndComboBox, CB_SETCURSEL, (WPARAM)2, (LPARAM)0);

  // Get aWindow of edit control in combobox created earlier.
  fHWndComboEditor = FindWindowEx(fHWndComboBox, nullptr, WC_EDIT, nullptr);

  //  Change the window procedure for the edit windows to the subclass procedure.
  origComboEditorWindowProc =
    (WNDPROC)SetWindowLongPtr(fHWndComboEditor, GWLP_WNDPROC, (LONG_PTR)ComboEditorWindowProc);

  // Create TreeView

  // Get the dimensions of the parent window's client area, and create
  // the tree-view control.
  GetClientRect(aWindow, &rcClient);
  fHWndHelpTree = CreateWindowEx(0, WC_TREEVIEW, TEXT("Tree View"),
    WS_VISIBLE | WS_CHILD | WS_BORDER | TVS_INFOTIP | TVS_HASBUTTONS | TVS_HASLINES |
      TVS_LINESATROOT,
    0, 0, rcClient.right, rcClient.bottom, aWindow, (HMENU)IDC_MAIN_TREE_VIEW,
    GetModuleHandle(nullptr), nullptr);

  // Initialize the Help Tree View.
  /*    if (!InitHelpTreeItems()) {
          DestroyWindow(fHWndHelpTree);
          return false;
      }*/

  // Create Status bar
  fHWndStatus = CreateWindowEx(0, STATUSCLASSNAME, nullptr, WS_CHILD | WS_VISIBLE | SBARS_SIZEGRIP,
    100, 100, 200, 200, aWindow, (HMENU)IDC_MAIN_STATUS, GetModuleHandle(nullptr), nullptr);

  SendMessage(fHWndStatus, SB_SETPARTS, sizeof(statwidths) / sizeof(int), (LPARAM)statwidths);
  // SendMessage(fHWndStatus, SB_SETTEXT, 0, (LPARAM) "Hi there :)");

  return true;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool G4UIWin32::ResizeComponents(HWND aWindow)
{
  RECT rcTool;
  G4int iToolHeight, iToolWidth;

  RECT rcStatus;
  G4int iStatusHeight;

  RECT rcComboBox;
  G4int iComboBoxHeight;

  G4int iTreeViewHeight, iTreeViewWidth;
  G4int iEditHeight, iEditWidth;

  RECT rcClient;

  // Size toolbar and get height and width
  SendMessage(fHWndToolBar, TB_AUTOSIZE, 0, 0);

  GetWindowRect(fHWndToolBar, &rcTool);
  iToolHeight = rcTool.bottom - rcTool.top;
  iToolWidth = rcTool.right - rcTool.left;

  // Size status bar and get height
  SendMessage(fHWndStatus, WM_SIZE, 0, 0);

  GetWindowRect(fHWndStatus, &rcStatus);
  iStatusHeight = rcStatus.bottom - rcStatus.top;

  // Size status the Combo Box and get height
  SendMessage(fHWndComboBox, WM_SIZE, 0, 0);

  GetWindowRect(fHWndComboBox, &rcComboBox);
  iComboBoxHeight = rcComboBox.bottom - rcComboBox.top;

  // Calculate remaining height and size edit
  GetClientRect(aWindow, &rcClient);

  iTreeViewHeight = rcClient.bottom - iToolHeight - iStatusHeight;
  iTreeViewWidth = iToolWidth / 4;

  iEditHeight = rcClient.bottom - iToolHeight - iComboBoxHeight - iStatusHeight;
  iEditWidth = iToolWidth - iTreeViewWidth;

  // TreeView location and size
  SetWindowPos(
    fHWndHelpTree, nullptr, 0, iToolHeight, iTreeViewWidth, iTreeViewHeight, SWP_NOZORDER);

  // Editor location and size
  SetWindowPos(
    fHWndEditor, nullptr, iTreeViewWidth, iToolHeight, iEditWidth, iEditHeight, SWP_NOZORDER);

  // ComboBox location and size
  SetWindowPos(
    fHWndComboBox, nullptr, iTreeViewWidth, iToolHeight + iEditHeight, iEditWidth, 200, 0);

  return true;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::ProcessTabKey()
{
  char buffer[256];

  // Only process the command it the user has written something
  if (SendMessage(fHWndComboBox, WM_GETTEXT, (WPARAM)sizeof(buffer), (LPARAM)buffer) != 0) {
    G4String command(buffer);

    SetFocus(fHWndComboBox);

    G4String cmd = Complete(command);
    const char* d = cmd.data();
    G4int l = strlen(d);
    Edit_SetText(fHWndComboEditor, d);
    Edit_SetSel(fHWndComboEditor, l, l);
  }
  else {
    if (GetFocus() == fHWndComboEditor)
      SetFocus(fHWndHelpTree);
    else
      SetFocus(fHWndComboBox);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::ProcessEscKey()
{
  // Clear the current selection.
  SendMessage(fHWndComboBox, CB_SETCURSEL, (WPARAM)(-1), (LPARAM)0);
  // Set the focus to the Combo Box.
  SetFocus(fHWndComboBox);
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::ProcessEnterKey()
{
  char buffer[256];
  DWORD dwIndex, numItems;

  // Only process the command it the user has written something
  if (SendMessage(fHWndComboBox, WM_GETTEXT, (WPARAM)sizeof(buffer), (LPARAM)buffer) != 0) {
    SetFocus(fHWndComboBox);

    // Read command
    G4String command(buffer);

    // Now clear the current selection.
    SendMessage(fHWndComboBox, CB_SETCURSEL, (WPARAM)-1, (LPARAM)0);

    if (fHelp) {
      exitHelp = true;
      fHelp = ConvertStringToInt(command.data(), fHelpChoice);
    }
    else {
      fHistory.push_back(command);
      fHistoryPos = -1;
      ApplyShellCommand(command, exitSession, exitPause);

      // Now update the history in the ComboBox

      // Check if this command exists in the ComboBox
      dwIndex = SendMessage(fHWndComboBox, CB_FINDSTRINGEXACT, (WPARAM)(-1), (LPARAM)buffer);

      //  Add the string, if necessary
      if (dwIndex == CB_ERR)
        dwIndex = SendMessage(fHWndComboBox, CB_INSERTSTRING, (WPARAM)0, (LPARAM)buffer);
      // If the string exists, move it to the first position
      if (dwIndex != CB_ERR) {
        SendMessage(fHWndComboBox, CB_DELETESTRING, (WPARAM)dwIndex, (LPARAM)0);
        dwIndex = SendMessage(fHWndComboBox, CB_INSERTSTRING, (WPARAM)0, (LPARAM)buffer);
      }

      numItems = SendMessage(fHWndComboBox, CB_GETCOUNT, (WPARAM)0, (LPARAM)0);
      while (numItems > MAX_HISTORY_ITEMS) {
        SendMessage(fHWndComboBox, CB_DELETESTRING, (WPARAM)(numItems - 1), (LPARAM)0);
        numItems = SendMessage(fHWndComboBox, CB_GETCOUNT, (WPARAM)0, (LPARAM)0);
      }
    }
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::ProcessUpKey()
{
  G4int pos = fHistoryPos == -1 ? fHistory.size() - 1 : fHistoryPos - 1;
  if ((pos >= 0) && (pos < (G4int)fHistory.size())) {
    G4String command = fHistory[pos];
    const char* d = command.data();
    G4int l = strlen(d);
    Edit_SetText(fHWndComboEditor, d);
    Edit_SetSel(fHWndComboEditor, l, l);

    fHistoryPos = pos;
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::ProcessDownKey()
{
  G4int pos = fHistoryPos + 1;
  if ((pos >= 0) && (pos < (G4int)fHistory.size())) {
    G4String command = fHistory[pos];
    const char* d = command.data();
    G4int l = strlen(d);
    Edit_SetText(fHWndComboEditor, d);
    Edit_SetSel(fHWndComboEditor, l, l);

    fHistoryPos = pos;
  }
  else if (pos >= (G4int)fHistory.size()) {
    Edit_SetText(fHWndComboEditor, "");
    Edit_SetSel(fHWndComboEditor, 0, 0);

    fHistoryPos = -1;
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool G4UIWin32::ProcessDefaultCommands(G4int idCommand)
{
  switch (idCommand) {
    case ID_EXIT_APP:
      PostMessage(fHWndMainWindow, WM_CLOSE, 0, 0);
      return true;
    case ID_OPEN_MACRO:
      DoOpenMacro(fHWndMainWindow);
      return true;
    case ID_SAVE_VIEWER_STATE:
      DoSaveViewer(fHWndMainWindow);
      return true;
    case ID_RUN_BEAMON:
      if (! fHelp) {
        G4String command = "/run/beamOn 1";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_RUN_CMD:
      return true;
    case ID_VIEW_SOLID:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/style s";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_VIEW_WIREFRAME:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/style w";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_PROJ_ORTHOGRAPHIC:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/projection o";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_PROJ_PERSPECTIVE:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/projection p";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_ZOOM_IN:
      if (! fHelp) {
        G4String command = "/vis/viewer/zoom 1.2";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_ZOOM_OUT:
      if (! fHelp) {
        G4String command = "/vis/viewer/zoom 0.8";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_ORIENTATION_XY:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/viewpointThetaPhi 0. 0.";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_ORIENTATION_XZ:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/viewpointThetaPhi 90. 0.";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_ORIENTATION_YZ:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/viewpointThetaPhi 0. 90.";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_ORIENTATION_OBLIQUE:
      if (! fHelp) {
        G4String command = "/vis/viewer/set/viewpointThetaPhi 45. -45.";
        ApplyShellCommand(command, exitSession, exitPause);
      }
      return true;
    case ID_HELP_ABOUT:
      return true;
    case ID_LOG_CLEAN:
      SetDlgItemText(fHWndMainWindow, IDC_MAIN_EDIT, "");
      return true;
    case ID_LOG_SAVE:
      DoSaveLog(fHWndMainWindow);
      return true;
    default:
      return false;
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4String G4UIWin32::GetToolTips(G4int idButton)
{
  switch (idButton) {
    case ID_OPEN_MACRO:
      return "Open and execute macro file";

    case ID_SAVE_VIEWER_STATE:
      return "Save viewer state";

    case ID_ZOOM_IN:
      return "Zoom in";

    case ID_ZOOM_OUT:
      return "Zoom out";

    case ID_RUN_BEAMON:
      return "Beam on (one particle)";

    case ID_HELP_ABOUT:
      return "About G4UIWin32";

    case ID_LOG_CLEAN:
      return "Clean log";

    case ID_LOG_SAVE:
      return "Save log";

    default:
      return "";
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4String G4UIWin32::GetHelpTreeToolTips(HTREEITEM item)
{
  // Tooltips for the help tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return "";
  G4UIcommandTree* treeTop = UI->GetTree();

  G4String itemText = GetItemPath(item);

  // Check if it is a command path
  if (TreeView_GetChild(fHWndHelpTree, item) != nullptr) itemText += "/";

  G4UIcommand* command = treeTop->FindPath(itemText.c_str());

  if (command) {
    // This is a command, return the first line of help
    return command->GetGuidanceLine(0).data();
  }
  else {
    // This is not a command, but a sub directory, return the title
    G4UIcommandTree* path = treeTop->FindCommandTree(itemText.c_str());
    if (path) return path->GetTitle().data();
  }

  return "";
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4String G4UIWin32::ConvertNewLines(G4String a_string)
{
  // Geant4 uses UNIX's style for new lines (\n)
  // we must convert them to Windows' style (\r\n)
  G4String str = std::move(a_string);
  size_t index = str.find("\n", 0);
  while (index < str.length()) {
    str.replace(index, 2, "\r\n");
    // Advance index forward so the next iteration doesn't pick it up as well.
    index = str.find("\n", index + 2);
  }
  return str;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::HelpTreeDoubleClick(HTREEITEM item)
{
  const char* item_path = GetItemPath(item);
  G4int l = strlen(item_path);
  Edit_SetText(fHWndComboEditor, item_path);
  Edit_SetSel(fHWndComboEditor, l, l);

  SetFocus(fHWndComboEditor);
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool G4UIWin32::SaveLogFile(LPCTSTR fileName)
{
  HANDLE hFile;
  G4bool bSuccess = false;

  hFile =
    CreateFile(fileName, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (hFile != INVALID_HANDLE_VALUE) {
    DWORD dwTextLength;

    dwTextLength = GetWindowTextLength(fHWndEditor);
    // No need to bother if there's no text.
    if (dwTextLength > 0) {
      LPSTR text;
      DWORD dwBufferSize = dwTextLength + 1;

      text = (LPSTR)GlobalAlloc(GPTR, dwBufferSize);
      if (text != nullptr) {
        if (GetWindowText(fHWndEditor, text, dwBufferSize)) {
          DWORD dwWritten;

          if (WriteFile(hFile, text, dwTextLength, &dwWritten, nullptr)) bSuccess = true;
        }
        GlobalFree(text);
      }
    }
    CloseHandle(hFile);
  }
  return bSuccess;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::AddText(LPSTR text)
{
  if ((fHWndEditor != nullptr) && (text != nullptr) && (text[0] != '\0')) {
    // Get current text length
    G4int ndx = GetWindowTextLength(fHWndEditor);

    // Select the end of the text
    SendMessage(fHWndEditor, EM_SETSEL, (WPARAM)ndx, (LPARAM)ndx);
    // Add the new text
    SendMessage(fHWndEditor, EM_REPLACESEL, (WPARAM)0, (LPARAM)text);
    // Scroll to the bottom
    SendMessage(fHWndEditor, WM_VSCROLL, SB_BOTTOM, NULL);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::DoOpenMacro(HWND aWindow)
{
  OPENFILENAME ofn;
  char szFileName[MAX_PATH] = "";

  ZeroMemory(&ofn, sizeof(ofn));

  ofn.lStructSize = sizeof(ofn);
  ofn.hwndOwner = aWindow;
  ofn.lpstrFilter = "Macro Files (*.mac)\0*.mac\0All Files (*.*)\0*.*\0";
  ofn.lpstrFile = szFileName;
  ofn.nMaxFile = MAX_PATH;
  ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
  ofn.lpstrDefExt = "mac";

  if (GetOpenFileName(&ofn)) {
    G4String command = "/control/execute " + G4String(szFileName);
    ApplyShellCommand(command, exitSession, exitPause);

    SendDlgItemMessage(aWindow, IDC_MAIN_STATUS, SB_SETTEXT, 0, (LPARAM) "Opened macro...");
    SendDlgItemMessage(aWindow, IDC_MAIN_STATUS, SB_SETTEXT, 1, (LPARAM)szFileName);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::DoSaveViewer(HWND aWindow)
{
  OPENFILENAME ofn;
  char szFileName[MAX_PATH] = "";

  ZeroMemory(&ofn, sizeof(ofn));

  ofn.lStructSize = sizeof(ofn);
  ofn.hwndOwner = aWindow;
  ofn.lpstrFilter = "Macro Files (*.mac)\0*.mac\0All Files (*.*)\0*.*\0";
  ofn.lpstrFile = szFileName;
  ofn.nMaxFile = MAX_PATH;
  ofn.lpstrDefExt = "mac";
  ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT;

  if (GetSaveFileName(&ofn)) {
    G4String command = "/vis/viewer/save " + G4String(szFileName);
    ApplyShellCommand(command, exitSession, exitPause);

    SendDlgItemMessage(aWindow, IDC_MAIN_STATUS, SB_SETTEXT, 0, (LPARAM) "State saved...");
    SendDlgItemMessage(aWindow, IDC_MAIN_STATUS, SB_SETTEXT, 1, (LPARAM)szFileName);
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::DoSaveLog(HWND aWindow)
{
  OPENFILENAME ofn;
  char szFileName[MAX_PATH] = "";

  ZeroMemory(&ofn, sizeof(ofn));

  ofn.lStructSize = sizeof(ofn);
  ofn.hwndOwner = aWindow;
  ofn.lpstrFilter = "Log Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0";
  ofn.lpstrFile = szFileName;
  ofn.nMaxFile = MAX_PATH;
  ofn.lpstrDefExt = "txt";
  ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT;

  if (GetSaveFileName(&ofn)) {
    if (SaveLogFile(szFileName)) {
      SendDlgItemMessage(aWindow, IDC_MAIN_STATUS, SB_SETTEXT, 0, (LPARAM) "Saved log file...");
      SendDlgItemMessage(aWindow, IDC_MAIN_STATUS, SB_SETTEXT, 1, (LPARAM)szFileName);
    }
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool G4UIWin32::InitHelpTreeItems()
{
  HTREEITEM newItem;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return false;
  G4UIcommandTree* treeTop = UI->GetTree();

  G4int treeSize = treeTop->GetTreeEntry();
  G4String commandText;
  for (G4int a = 0; a < treeSize; a++) {
    // Creating new item
    commandText = treeTop->GetTree(a + 1)->GetPathName().data();

    // Add the item to the tree-view control.
    newItem = AddItemToHelpTree(const_cast<LPTSTR>(GetShortCommandPath(commandText).c_str()));

    if (newItem == nullptr) return false;

    // Look for children
    CreateHelpTree(newItem, treeTop->GetTree(a + 1));
  }

  return true;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4UIWin32::CreateHelpTree(HTREEITEM aParent, G4UIcommandTree* aCommandTree)
{
  if ((aParent != nullptr) && (aCommandTree != nullptr)) {
    // Creating new item
    HTREEITEM newItem;

    G4String commandText;
    // Get the Sub directories
    for (G4int a = 0; a < aCommandTree->GetTreeEntry(); a++) {
      commandText = aCommandTree->GetTree(a + 1)->GetPathName().data();

      // Add the item to the tree-view control.
      newItem =
        AddItemToHelpTree(const_cast<LPTSTR>(GetShortCommandPath(commandText).c_str()), aParent);

      // Look for children
      CreateHelpTree(newItem, aCommandTree->GetTree(a + 1));
    }

    // Get the Commands
    for (G4int a = 0; a < aCommandTree->GetCommandEntry(); a++) {
      commandText = aCommandTree->GetCommand(a + 1)->GetCommandPath().data();

      // Add the item to the tree-view control.
      AddItemToHelpTree(const_cast<LPTSTR>(GetShortCommandPath(commandText).c_str()), aParent);
    }
  }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
HTREEITEM G4UIWin32::AddItemToHelpTree(LPTSTR lpszItem, HTREEITEM aParent)
{
  TVITEM tvi;
  TVINSERTSTRUCT tvins;
  static auto hPrev = (HTREEITEM)TVI_FIRST;

  tvi.mask = TVIF_TEXT | TVIF_IMAGE | TVIF_SELECTEDIMAGE | TVIF_PARAM;

  // Set the text of the item.
  tvi.pszText = lpszItem;
  tvi.cchTextMax = sizeof(tvi.pszText) / sizeof(tvi.pszText[0]);

  // Save the heading level in the item's application-defined
  // data area.
  tvi.lParam = (LPARAM)aParent;
  tvins.item = tvi;
  tvins.hInsertAfter = hPrev;
  // Set the parent item.
  tvins.hParent = aParent;

  // Add the item to the tree-view control.
  hPrev = (HTREEITEM)SendMessage(
    fHWndHelpTree, TVM_INSERTITEM, (WPARAM)0, (LPARAM)(LPTVINSERTSTRUCT)&tvins);

  if (hPrev == nullptr)
    return nullptr;
  else
    return hPrev;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4String G4UIWin32::GetShortCommandPath(G4String commandPath)
{
  G4String str = std::move(commandPath);

  if (str.find_last_of("/") == (str.size() - 1)) str = str.erase(str.size() - 1, 1);

  str = str.erase(0, str.find_last_of("/") + 1);

  return str;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
LPSTR G4UIWin32::GetItemPath(HTREEITEM item)
{
  // Get the text for the item.
  TVITEM tvitem;
  tvitem.mask = TVIF_TEXT;
  tvitem.hItem = item;
  TCHAR infoTipBuf[1024];
  tvitem.pszText = infoTipBuf;
  tvitem.cchTextMax = sizeof(infoTipBuf) / sizeof(TCHAR);

  std::string str = "";
  while (item != nullptr) {
    TreeView_GetItem(fHWndHelpTree, &tvitem);
    str = "/" + std::string(tvitem.pszText) + str;

    item = TreeView_GetParent(fHWndHelpTree, item);
    tvitem.hItem = item;
  }

  auto* result = new TCHAR[str.size() + 1];
  result[str.size()] = 0;
  std::copy(str.begin(), str.end(), result);

  return result;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/

/****************************************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
G4bool ConvertStringToInt(const char* aString, G4int& aInt)
{
  aInt = 0;
  if (aString == nullptr) return false;
  char* s;
  G4long value = strtol(aString, &s, 10);
  if (s == aString) return false;
  aInt = value;
  return true;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/****************************************************************************************************/
