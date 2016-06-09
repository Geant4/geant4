//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Win32.cc,v 1.8 2004/04/08 15:35:05 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// G.Barrand

#if defined(G4INTY_BUILD_WIN32) || defined(G4INTY_USE_WIN32)

// this :
#include "G4Win32.hh"

#include "G4ios.hh"

static char className[] = "G4Win32";

G4Win32* G4Win32::instance  = NULL;

static G4bool    Win32Inited   = FALSE;
static HWND      topWindow     = NULL;
/***************************************************************************/
G4Win32* G4Win32::getInstance (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if (instance==NULL) {
    instance = new G4Win32();
  }
  return instance;
}
/***************************************************************************/
G4Win32::G4Win32 (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(Win32Inited==FALSE) { // Should be Done once.

    WNDCLASS         wc;
    wc.style         = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc   = (WNDPROC)DefWindowProc;
    wc.cbClsExtra    = 0;
    wc.cbWndExtra    = 0;
    wc.hInstance     = ::GetModuleHandle(NULL);
    wc.hIcon         = LoadIcon  (NULL,IDI_APPLICATION);
    wc.hCursor       = LoadCursor(NULL,IDC_ARROW);
    wc.hbrBackground = GetStockBrush(BLACK_BRUSH);
    wc.lpszMenuName  = className;
    wc.lpszClassName = className;
    ::RegisterClass  (&wc);
    
    topWindow   = ::CreateWindow(className,className, 
				 WS_OVERLAPPEDWINDOW,
				 CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, 
				 NULL, NULL, 
				 ::GetModuleHandle(NULL),
				 NULL);
    
    if(topWindow==NULL) {
      G4cout << "G4Win32 : Unable to create Win32 window." << G4endl;
    }

    Win32Inited = TRUE;
  }

  AddDispatcher((G4DispatchFunction)G4Win32::dispatchWin32Event);
  SetMainInteractor(topWindow);
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
  MSG event;
  while ( ::PeekMessage(&event, NULL, 0, 0, PM_REMOVE) ) {
    ::TranslateMessage(&event);
    ::DispatchMessage (&event);
  }
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

