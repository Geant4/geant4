// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Win32.cc,v 1.3 1999-05-07 10:34:18 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

#if defined(G4INTY_BUILD_WIN32) || defined(G4INTY_USE_WIN32)

#include <stdlib.h>
#include "G4ios.hh"

#include "G4UImanager.hh"

#include "G4Win32.hh"

static char className[] = "G4Win32";

G4Win32* G4Win32::instance  = NULL;

static G4bool    Win32Inited   = FALSE;
static HINSTANCE hInstance     = NULL;
static HINSTANCE hPrevInstance = NULL;
static LPTSTR    lpszCmdLine   = NULL;
static int       nCmdShow      = 0;
static HWND      topWindow     = NULL;
/***************************************************************************/
G4Win32* G4Win32::getInstance (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return G4Win32::getInstance(NULL,NULL,NULL,0);
}
/***************************************************************************/
G4Win32* G4Win32::getInstance (
 HINSTANCE a_hInstance
,HINSTANCE a_hPrevInstance
,LPSTR  a_lpszCmdLine
,int    a_nCmdShow
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if (instance==NULL) {
    instance = new G4Win32(a_hInstance,
			   a_hPrevInstance,
			   a_lpszCmdLine,
			   a_nCmdShow);
  }
  return instance;
}
/***************************************************************************/
G4Win32::G4Win32 (
 HINSTANCE a_hInstance
,HINSTANCE a_hPrevInstance
,LPSTR  a_lpszCmdLine
,int    a_nCmdShow
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(Win32Inited==FALSE) { // Should be Done once.

    hInstance        = a_hInstance;
    hPrevInstance    = a_hPrevInstance;
    lpszCmdLine      = a_lpszCmdLine;
    nCmdShow         = a_nCmdShow;
  
    if(hInstance==NULL) {
      G4cout << "G4Win32::G4Win32 : NULL hInstance given." <<endl;
    }

    if(hPrevInstance==NULL) {
      WNDCLASS         wc;
      wc.style         = CS_HREDRAW | CS_VREDRAW;
      wc.lpfnWndProc   = (WNDPROC)DefWindowProc;
      wc.cbClsExtra    = 0;
      wc.cbWndExtra    = 0;
      wc.hInstance     = hInstance;
      wc.hIcon         = LoadIcon  (NULL,IDI_APPLICATION);
      wc.hCursor       = LoadCursor(NULL,IDC_ARROW);
      wc.hbrBackground = GetStockBrush(BLACK_BRUSH);
      wc.lpszMenuName  = className;
      wc.lpszClassName = className;
      ::RegisterClass  (&wc);
    }
    
    topWindow   = ::CreateWindow(className,className, 
				 WS_OVERLAPPEDWINDOW,
				 CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, 
				 NULL, NULL, hInstance, NULL);
    
    if(topWindow==NULL) {
      G4cout << "G4Win32 : Unable to create Win32 window." << endl;
    }

    Win32Inited = TRUE;
  }

  AddDispatcher     ((G4DispatchFunction)G4Win32::dispatchWin32Event);
  SetMainInteractor (topWindow);
}
/***************************************************************************/
G4Win32::~G4Win32 (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(this==instance) {
    instance = NULL;
  }
}
/***************************************************************************/
G4bool G4Win32::Inited (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return Win32Inited;
}
/***************************************************************************/
void* G4Win32::GetEvent (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  static MSG event;
  BOOL status = ::GetMessage(&event, NULL, 0, 0);
  if(status==FALSE) return NULL;
  return &event;
}
/***************************************************************************/
void G4Win32::FlushAndWaitExecution (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4Win32::getWinMainArguments (
 HINSTANCE* a_hInstance
,HINSTANCE* a_hPrevInstance
,LPSTR*  a_lpszCmdLine
,int*    a_nCmdShow
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
 if(a_hInstance!=NULL)     *a_hInstance     = hInstance;
 if(a_hPrevInstance!=NULL) *a_hPrevInstance = hPrevInstance;
 if(a_lpszCmdLine!=NULL)   *a_lpszCmdLine   = lpszCmdLine;
 if(a_nCmdShow!=NULL)      *a_nCmdShow      = nCmdShow;
}
/***************************************************************************/
G4bool G4Win32::dispatchWin32Event  (
 void* a_event
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  ::TranslateMessage((MSG*)a_event);
  ::DispatchMessage ((MSG*)a_event);
  return TRUE;
}

#endif //HAS_WIN32

