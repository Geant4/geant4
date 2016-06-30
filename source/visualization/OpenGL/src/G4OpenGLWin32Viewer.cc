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
// $Id: G4OpenGLWin32Viewer.cc 97430 2016-06-03 07:20:11Z gcosmo $
//
// 
// G4OpenGLWin32Viewer : Class to provide WindowsNT specific
//                     functionality for OpenGL in GEANT4
//
// 27/06/2003 : G.Barrand : implementation (at last !).

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLWin32Viewer.hh"
#include "G4VViewer.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"

#include "G4ios.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"


//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::SetView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fHDC) return;
  if(!fHGLRC) return;
  ::wglMakeCurrent(fHDC,fHGLRC);
  G4OpenGLViewer::SetView ();  
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::ShowView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fHDC) return;
  glFlush ();
  // Empty the Windows message queue :
  MSG event;
  while ( ::PeekMessage(&event, NULL, 0, 0, PM_REMOVE) ) {
    ::TranslateMessage(&event);
    ::DispatchMessage (&event);
  }
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::GetWin32Connection (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::CreateGLWin32Context (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::CreateMainWindow (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWindow) return; //Done.

  // Bill Gates stuff...
  static const char className[] = "G4OpenGLWin32";
  static bool done = false;
  if(done==false) {
    WNDCLASS wc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = (WNDPROC)WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = ::GetModuleHandle(NULL);
    wc.hIcon = LoadIcon  (NULL, IDI_APPLICATION);
    wc.hCursor = LoadCursor(NULL,IDC_CROSS);
    wc.hbrBackground = NULL;
    wc.lpszMenuName = className;
    wc.lpszClassName = className;
    ::RegisterClass(&wc);
    done = true;
  }  
  
  ResizeWindow(fVP.GetWindowSizeHintX(),fVP.GetWindowSizeHintY());

  int x_res=GetSystemMetrics(SM_CXSCREEN);
  int y_res=GetSystemMetrics(SM_CYSCREEN);
  
  //FIXME : NOT tested !
  fWindow = ::CreateWindow(className,fName.c_str(), 
			   WS_OVERLAPPEDWINDOW,
			   //WS_CHILD | WS_VISIBLE,
                           //			   0,0,
                           fVP.GetWindowAbsoluteLocationHintX(x_res),
                           fVP.GetWindowAbsoluteLocationHintY(y_res),
			   getWinWidth(), getWinHeight(),
			   NULL, NULL, 
			   ::GetModuleHandle(NULL),
			   NULL);
  if(!fWindow) return;

  ::SetWindowLongPtr(fWindow,GWLP_USERDATA,LONG_PTR(this));

  // initialize OpenGL rendering :
  fHDC = ::GetDC(fWindow);
  if( fHDC && (SetWindowPixelFormat(fHDC)==TRUE) ) {
    fHGLRC = ::wglCreateContext(fHDC);
  }
  
  if(fHDC && fHGLRC) {
    ::wglMakeCurrent(fHDC,fHGLRC);
  }

  ::SetForegroundWindow(fWindow);
  ::ShowWindow(fWindow,SW_SHOWDEFAULT);
  ::UpdateWindow(fWindow);
  ::DrawMenuBar(fWindow);
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLWin32Viewer::G4OpenGLWin32Viewer (
 G4OpenGLSceneHandler& scene
)
:G4VViewer (scene, -1)
,G4OpenGLViewer (scene)
,fWindow(0)
,fHDC(0)
,fHGLRC(0)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLWin32Viewer::~G4OpenGLWin32Viewer (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  // This is the end (Jim Morisson).
  if (fViewId >= 0) {
    if(wglGetCurrentContext()!=NULL) wglMakeCurrent(NULL,NULL);
    if(fHGLRC)	{
      wglDeleteContext(fHGLRC);
      fHGLRC = NULL;
    }
    
    if(fWindow) {
      ::SetWindowLongPtr(fWindow,GWLP_USERDATA,LONG(NULL));
      if(fHDC) ::ReleaseDC(fWindow,fHDC);
      ::DestroyWindow(fWindow);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
LRESULT CALLBACK G4OpenGLWin32Viewer::WindowProc ( 
 HWND   aWindow
,UINT   aMessage
,WPARAM aWParam
,LPARAM aLParam
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
/*
  switch (aMessage) { 
  case WM_PAINT:{
    printf("debug : PAINT\n");
    HDC	hDC;
    PAINTSTRUCT	ps;
    hDC = BeginPaint(aWindow,&ps);
    if(This) {
      // FIXME : To have an automatic refresh someone have to redraw here.
    }
    EndPaint(aWindow, &ps);

    //FIXME : have to handle WM_RESIZE
    //pView->fWinSize_x = (G4int) width;
    //pView->fWinSize_y = (G4int) height;
    G4OpenGLWin32Viewer* This = 
      (G4OpenGLWin32Viewer*)::GetWindowLong(aWindow,GWL_USERDATA);
    if(This) {
      This->SetView();
      glViewport(0,0,This->fWinSize_x,This->fWinSize_y);
      This->ClearView();
      This->DrawView();
      // WARNING : the below empty the Windows message queue...
      This->FinishView();
    }
  } return 0;
  default:
    return DefWindowProc(aWindow,aMessage,aWParam,aLParam);
  }
*/
  return DefWindowProc(aWindow,aMessage,aWParam,aLParam);
}

//////////////////////////////////////////////////////////////////////////////
bool G4OpenGLWin32Viewer::SetWindowPixelFormat(
 HDC aHdc
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  // The ungessable...

  PIXELFORMATDESCRIPTOR pfd;
  pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);
  pfd.nVersion = 1;
  pfd.dwFlags = 
    PFD_DRAW_TO_WINDOW | 
    PFD_SUPPORT_OPENGL | 
    PFD_DOUBLEBUFFER | 
    PFD_STEREO_DONTCARE;  
  pfd.iPixelType = PFD_TYPE_RGBA;
  pfd.cColorBits = 32;
  pfd.cRedBits = 8;
  pfd.cRedShift = 16;
  pfd.cGreenBits = 8;
  pfd.cGreenShift = 8;
  pfd.cBlueBits	= 8;
  pfd.cBlueShift = 0;
  pfd.cAlphaBits = 0;
  pfd.cAlphaShift = 0;
  pfd.cAccumBits = 64;	
  pfd.cAccumRedBits = 16;
  pfd.cAccumGreenBits = 16;
  pfd.cAccumBlueBits = 16;
  pfd.cAccumAlphaBits = 0;
  pfd.cDepthBits = 32;
  pfd.cStencilBits = 8;
  pfd.cAuxBuffers = 0;
  pfd.iLayerType = PFD_MAIN_PLANE;
  pfd.bReserved	= 0;
  pfd.dwLayerMask = 0;
  pfd.dwVisibleMask = 0;
  pfd.dwDamageMask = 0;
  
  int pixelIndex = ::ChoosePixelFormat(aHdc,&pfd);
  if (pixelIndex==0) {
    pixelIndex = 1;	
    if (::DescribePixelFormat(aHdc, 
			      pixelIndex, 
			      sizeof(PIXELFORMATDESCRIPTOR), 
			      &pfd)==0) {
      return false;
    }
  }
  if (::SetPixelFormat(aHdc,pixelIndex,&pfd)==FALSE) return false;
  return true;
}


#endif
