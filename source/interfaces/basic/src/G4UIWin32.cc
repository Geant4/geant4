// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIWin32.cc,v 1.3 1999-04-16 10:06:00 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

//#define DEBUG

#ifdef G4UI_BUILD_WIN32_SESSION

#include <string.h>

#include <windows.h>
#include <windowsx.h>
#include <wingdi.h>

#include <Strstrea.h>

#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommandStatus.hh"
#include "G4Win32.hh"

#include "G4UIWin32.hh"

#define TEXT_MAX_LINES 300

class TextBuffer {
public:
  TextBuffer():linei(0),linen(TEXT_MAX_LINES),endOfPage(0),heightOfPage(12) {
    lines = new G4String[linen];
  }
  ~TextBuffer() {
    delete [] lines;
  }
  int GetNumberOfLines () {
    return linei;
  }
  void SetHeightOfPage (int a_height) {
    heightOfPage = a_height;
  }
  void SetEndOfPage (int a_value) {
    if( (a_value<0) || (a_value>=linei)) {
      endOfPage = linei-1;
    } else {
      endOfPage = a_value;
    }
  }
  int GetEndOfPage () {
    return endOfPage;
  }
  void IncrementEndOfPage () {
    endOfPage++;
    if(endOfPage>=linei) endOfPage = linei-1;
  }
  void DecrementEndOfPage () {
    endOfPage--;
    if(endOfPage<0) endOfPage = 0;
  }
  void JumpDownEndOfPage () {
    endOfPage += heightOfPage;
    if(endOfPage>=linei) endOfPage = linei-1;
  }
  void JumpUpEndOfPage () {
    endOfPage -= heightOfPage;
    if(endOfPage<0) endOfPage = 0;
  }
  G4bool AppendString (char* a_string) {
    G4bool value = false;
    if( (a_string==NULL) || (a_string[0]=='\0') ) return value;
    int length = strlen(a_string);
    if(a_string[length-1]=='\n') {
      lines[linei] += a_string; 
      lines[linei] = lines[linei].strip(G4String::trailing,'\n');
      linei++;
      value = true;
    } else {
      lines[linei] += a_string; 
    }
    if(linei>=linen) {
      for(int count=0;count<linen;count++) {
	lines[count] = "";
      }
      linei = 0;
    }
    if(value==true) endOfPage = linei-1;
    return value;
  }
  void Draw (HDC a_hdc,RECT* a_rect) {
    TEXTMETRIC tm;
    GetTextMetrics (a_hdc,&tm);
    short charWidth = tm.tmAveCharWidth;
    short charHeight = tm.tmHeight + tm.tmExternalLeading;
    for(int row=0;row<heightOfPage;row++) {
      int rowi = endOfPage - row;
      short y = a_rect->bottom - charHeight * (row + 1);
      if((rowi>=0)&&(rowi<linei)) {
	const char* string = lines[rowi].data();
	if(string!=NULL) {
	  TextOut (a_hdc,0,y,(char*)string,strlen((char*)string));
	}
      }
    }
  }
private:
  G4String* lines;
  int linen;
  int linei;
  int endOfPage,heightOfPage;
};

static char mainClassName[] = "G4UIWin32";
static char textClassName[] = "G4UIWin32/Text";
static G4bool exitSession = true;
static G4bool exitPause = true;
static G4bool exitHelp = true;
static G4UIsession* tmpSession = NULL;

static FARPROC oldEditWindowProc;
static G4bool ConvertStringToInt(const char*,int&);

static int actionIdentifier = 0;

static unsigned CommandsHashFun(const int& identifier) {
  return (unsigned)identifier;
}
/***************************************************************************/
G4UIWin32::G4UIWin32 (
 HANDLE a_hInstance
,HANDLE a_hPrevInstance
,LPSTR  a_lpszCmdLine
,int    a_nCmdShow
)
:mainWindow(NULL)
,textWindow(NULL)
,editWindow(NULL)
,menuBar(NULL)
,commands(CommandsHashFun)
,textBuffer(NULL)
,textCols(80)
,textRows(12)
,fHelp(false)
,fHelpChoice(0)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  interactorManager = G4Win32::getInstance (a_hInstance,a_hPrevInstance,
					    a_lpszCmdLine,a_nCmdShow);
  static G4bool Done = FALSE;
  if(Done==FALSE) {
    if(a_hPrevInstance==NULL) {
      WNDCLASS         wc;
      wc.style         = CS_HREDRAW | CS_VREDRAW;
      wc.lpfnWndProc   = (WNDPROC)G4UIWin32::MainWindowProc;
      wc.cbClsExtra    = 0;
      wc.cbWndExtra    = 0;
      wc.hInstance     = a_hInstance;
      wc.hIcon         = LoadIcon(NULL,IDI_APPLICATION);
      wc.hCursor       = LoadCursor(NULL,IDC_ARROW);
      wc.hbrBackground = GetStockObject(BLACK_BRUSH);
      wc.lpszMenuName  = mainClassName;
      wc.lpszClassName = mainClassName;
      ::RegisterClass  (&wc);

      wc.style         = CS_HREDRAW | CS_VREDRAW;
      wc.lpfnWndProc   = (WNDPROC)G4UIWin32::TextWindowProc;
      wc.cbClsExtra    = 0;
      wc.cbWndExtra    = 0;
      wc.hInstance     = a_hInstance;
      wc.hIcon         = LoadIcon(NULL,IDI_APPLICATION);
      wc.hCursor       = LoadCursor(NULL,IDC_ARROW);
      wc.hbrBackground = GetStockObject(WHITE_BRUSH);
      wc.lpszMenuName  = textClassName;
      wc.lpszClassName = textClassName;
      ::RegisterClass  (&wc);
    }
    Done = TRUE;
  }

  /*
  */
  menuBar = CreateMenu();
  defaultMenu = CreatePopupMenu();
  AppendMenu(menuBar,MF_POPUP,(UINT)defaultMenu,"Geant4");

  textBuffer = new TextBuffer();

  tmpSession = this;
  mainWindow = ::CreateWindow(mainClassName,mainClassName, 
			   WS_OVERLAPPEDWINDOW,
			   CW_USEDEFAULT,CW_USEDEFAULT, 
			   0,0,
			   NULL,menuBar,a_hInstance,NULL);
  tmpSession = NULL;
  ::SetWindowLong(mainWindow,GWL_USERDATA,LONG(this));

  ::SetForegroundWindow(mainWindow);
  ::ShowWindow(mainWindow,SW_SHOWDEFAULT);
  ::UpdateWindow(mainWindow);

  if(UI!=NULL) UI->SetCoutDestination(this);
}
/***************************************************************************/
G4UIWin32::~G4UIWin32 (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{ 
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) {
    UI->SetSession(NULL);
    UI->SetCoutDestination(NULL);
  }
  delete textBuffer;
  ::DestroyWindow(mainWindow);
}
/***************************************************************************/
G4UIsession* G4UIWin32::SessionStart (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(interactorManager==NULL) return this;
  Prompt       ("session");
  exitSession  = false;
  interactorManager->DisableSecondaryLoop ();
  void*         event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitSession==true) break;
  }
  interactorManager->EnableSecondaryLoop ();
  return       this;
}
/***************************************************************************/
void G4UIWin32::Prompt (
 G4String a_prompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4UIWin32::SessionTerminate (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4UIWin32::PauseSessionStart (
 G4String a_state
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_state=="G4_pause> ") { 
    SecondaryLoop ("Pause, type continue to exit this state");
  }

  if(a_state=="EndOfEvent") {
    // Picking with feed back in event data Done here !!!
    SecondaryLoop ("End of event, type continue to exit this state");
  }
}
/***************************************************************************/
void G4UIWin32::SecondaryLoop (
 G4String a_prompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(interactorManager==NULL) return;
  Prompt(a_prompt);
  exitPause = false;
  void* event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitPause==true) break;
  }
  Prompt("session");
}
/***************************************************************************/
G4int G4UIWin32::ReceiveG4cout (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  TextAppendString((char*)a_string.data());
  return 0;
}
/***************************************************************************/
G4int G4UIWin32::ReceiveG4cerr (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  TextAppendString((char*)a_string.data());
  return 0;
}
/***************************************************************************/
G4bool G4UIWin32::GetHelpChoice(
 G4int& aInt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  fHelp = true;
  //
  if(interactorManager==NULL) return false;
  Prompt("Help");
  exitHelp = false;
  void* event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitHelp==true) break;
  }
  Prompt("session");
  //
  if(fHelp==false) return false;
  aInt = fHelpChoice;
  fHelp = false;
  return true;
}
/***************************************************************************/
void G4UIWin32::ExitHelp(
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4UIWin32::AddMenu (
 const char* a_name
,const char* a_label
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_name==NULL) return;
  if(defaultMenu!=NULL) {
    DeleteMenu (menuBar,0,MF_BYPOSITION);
    defaultMenu = NULL;
  }
  HMENU hMenu = CreatePopupMenu();
  AppendMenu(menuBar,MF_POPUP,(UINT)hMenu,a_label);
  AddInteractor(a_name,(G4Interactor)hMenu);
  DrawMenuBar(mainWindow);
}
/***************************************************************************/
void G4UIWin32::AddButton (
 const char* a_menu
,const char* a_label
,const char* a_command
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_menu==NULL) return;
  if(a_label==NULL) return;
  if(a_command==NULL) return;
  HMENU hMenu = (HMENU)GetInteractor(a_menu);
  actionIdentifier++;
  commands[actionIdentifier] = a_command;
  AppendMenu (hMenu,MF_STRING,actionIdentifier,a_label);
}
/***************************************************************************/
G4String G4UIWin32::GetCommand (
 int a_id
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return commands[a_id];
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
LRESULT CALLBACK G4UIWin32::MainWindowProc ( 
 HWND   a_window
,UINT   a_message
,WPARAM a_wParam
,LPARAM a_lParam
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  static short charWidth,charHeight;

  switch (a_message) { 
  case WM_CREATE:{
    HDC hdc;
    TEXTMETRIC tm;
    RECT rect;
    GetWindowRect (a_window,&rect);

    hdc = GetDC (a_window);
    GetTextMetrics (hdc,&tm);
    charWidth = tm.tmAveCharWidth;
    charHeight = tm.tmHeight + tm.tmExternalLeading;
    ReleaseDC (a_window,hdc);

    G4UIWin32* This = (G4UIWin32*)tmpSession;

    This->textWindow = CreateWindow (textClassName,NULL,
			 WS_CHILD | WS_VISIBLE | WS_VSCROLL,
			 0,0,
			 This->textCols * charWidth,
			 This->textRows * charHeight,
			 a_window,NULL,
			 GetWindowInstance(a_window),
			 NULL);
    ::SetWindowLong (This->textWindow,GWL_USERDATA,LONG(This));

    This->editWindow = CreateWindow ("edit",NULL,
			 WS_CHILD | WS_VISIBLE | WS_BORDER,
			 0,This->textRows  * charHeight,
			 This->textCols  * charWidth,charHeight,
			 a_window,(HMENU)1,
			 GetWindowInstance(a_window),
			 NULL);
    oldEditWindowProc = (FARPROC)GetWindowLong(This->editWindow,GWL_WNDPROC);
    SetWindowLong (This->editWindow,GWL_WNDPROC,(LONG)EditWindowProc);

    MoveWindow (a_window,
		rect.left,rect.top,
		2 * GetSystemMetrics(SM_CXFRAME) + 
		This->textCols  * charWidth,
		GetSystemMetrics(SM_CYCAPTION) + 
		2 * GetSystemMetrics(SM_CYFRAME) + 
		This->textRows * charHeight + charHeight,
		TRUE);

    }return 0;
  case WM_SIZE:{
    G4UIWin32* This = (G4UIWin32*)::GetWindowLong(a_window,GWL_USERDATA);
    if(This!=NULL) {
      // Client size :
      int width = LOWORD(a_lParam);
      int height = HIWORD(a_lParam);
      int editHeight = /*2 * GetSystemMetrics(SM_CYBORDER) + */ charHeight;
      MoveWindow (This->textWindow,
		  0,0,
		  width,height - editHeight,
		  FALSE);
      MoveWindow (This->editWindow,
		  0,height - editHeight,
		  width,charHeight,
		  FALSE);
      ((TextBuffer*)This->textBuffer)->SetHeightOfPage(height/charHeight);
    }
    }return 0;
  case WM_SETFOCUS:{
    G4UIWin32* This = (G4UIWin32*)::GetWindowLong(a_window,GWL_USERDATA);
    SetFocus (This->editWindow);
    }return 0;
  case WM_COMMAND:{
    G4UIWin32* This = (G4UIWin32*)::GetWindowLong(a_window,GWL_USERDATA);
    if(This!=NULL) {
      if(This->fHelp==false) {
	G4String command = This->GetCommand(a_wParam);
	This->ApplyShellCommand (command,exitSession,exitPause);
      }
    }
    }return 0;
  case WM_DESTROY:
    PostQuitMessage(0);
    return   0;
  }
  return (DefWindowProc(a_window,a_message,a_wParam,a_lParam));
}
/***************************************************************************/
LRESULT CALLBACK G4UIWin32::TextWindowProc ( 
 HWND   a_window
,UINT   a_message
,WPARAM a_wParam
,LPARAM a_lParam
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  switch (a_message) { 
  case WM_PAINT:{
    G4UIWin32* This = (G4UIWin32*)::GetWindowLong(a_window,GWL_USERDATA);
    if(This!=NULL) {
      TextBuffer* textBuffer = (TextBuffer*)This->textBuffer;
      RECT rect;
      GetClientRect (a_window,&rect);
      PAINTSTRUCT ps;
      HDC hdc = BeginPaint(a_window,&ps);
      textBuffer->Draw(hdc,&rect);
      EndPaint(a_window,&ps);
    }
    }return 0;
  case WM_VSCROLL:{
    G4UIWin32* This = (G4UIWin32*)::GetWindowLong(a_window,GWL_USERDATA);
    if(This!=NULL) {
      TextBuffer* textBuffer = (TextBuffer*)This->textBuffer;
      int what = LOWORD(a_wParam);
      switch(what) {
      case SB_LINEUP:
	textBuffer->DecrementEndOfPage();
	break;
      case SB_LINEDOWN:
	textBuffer->IncrementEndOfPage();
	break;
      case SB_PAGEUP:
	textBuffer->JumpUpEndOfPage();
	break;
      case SB_PAGEDOWN:
	textBuffer->JumpDownEndOfPage();
	break;
      case SB_THUMBPOSITION:
      case SB_THUMBTRACK:
	textBuffer->SetEndOfPage(HIWORD(a_wParam));
	break;
      default:
	return 0;
      }
      int eop = textBuffer->GetEndOfPage();
      SetScrollPos(a_window,SB_VERT,eop,TRUE);
      InvalidateRect(a_window,NULL,TRUE);
    }}return 0;
  case WM_DESTROY:
    PostQuitMessage(0);
    return 0;
  }
  return (DefWindowProc(a_window,a_message,a_wParam,a_lParam));
}
/***************************************************************************/
LRESULT CALLBACK G4UIWin32::EditWindowProc ( 
 HWND   a_window
,UINT   a_message
,WPARAM a_wParam
,LPARAM a_lParam
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  switch (a_message) { 
  case WM_KEYDOWN:
    switch(a_wParam){
    case VK_RETURN:{
      G4UIWin32* This = (G4UIWin32*)::GetWindowLong(
			 GetParent(a_window),GWL_USERDATA);
      char buffer[128];
      GetWindowText (a_window,buffer,128);
      G4String command (buffer);
      SetWindowText (a_window,"");

      if(This->fHelp==true) {
	exitHelp = true;
	This->fHelp = ConvertStringToInt(command.data(),This->fHelpChoice);
      } else {
	This->ApplyShellCommand (command,exitSession,exitPause);
      }

      }break;
    case VK_TAB:{
      G4UIWin32* This = (G4UIWin32*)::GetWindowLong(
			 GetParent(a_window),GWL_USERDATA);
      if(This->fHelp==true) break;
      char buffer[128];
      Edit_GetText(a_window,buffer,128);

      G4String command(buffer);

      const char* d = This->Complete(command).data();
      Edit_SetText(a_window,d);
      int pos = strlen(d);
      Edit_SetSel(a_window,pos,pos);

      }break;
    }
  }
  return CallWindowProc(oldEditWindowProc,
			a_window,a_message,
			a_wParam,a_lParam);
}
/***************************************************************************/
void G4UIWin32::TextAppendString (
 char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if( (a_string==NULL) || (a_string[0]=='\0') ) return;
  if(textWindow==NULL) return;
  if(((TextBuffer*)textBuffer)->AppendString(a_string)==true) {
    // The appending triggers and end of line, and then updates window :
    RECT rect;
    GetClientRect(textWindow,&rect);
    InvalidateRect(textWindow,NULL,TRUE); //To erase background.
    PAINTSTRUCT ps;
    HDC hdc = GetDC(textWindow);
    ((TextBuffer*)textBuffer)->Draw(hdc,&rect);
    ReleaseDC (textWindow,hdc);
    int linen = ((TextBuffer*)textBuffer)->GetNumberOfLines();
    SetScrollRange(textWindow,SB_VERT,0,linen-1,TRUE);
    SetScrollPos(textWindow,SB_VERT,linen-1,TRUE);
  }
}
//////////////////////////////////////////////////////////////////////////////
G4bool ConvertStringToInt(
 const char* aString
,int& aInt
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  aInt = 0;
  if(aString==NULL) return false;
  char* s;
  long value = strtol(aString,&s,10);
  if(s==aString) return false;
  aInt = value;
  return true;
}


#endif
