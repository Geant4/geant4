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
// G4OpenGLWin32Viewer : Class to provide Windows specific
//                       functionality for OpenGL
//
// 27/06/2003 : G.Barrand : implementation (at last !).

#include "G4OpenGLWin32Viewer.hh"
#include "G4VViewer.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4Scene.hh"

#include "G4ios.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

#include "G4SystemOfUnits.hh"

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
  static G4bool done = false;
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

  G4int x_res=GetSystemMetrics(SM_CXSCREEN);
  G4int y_res=GetSystemMetrics(SM_CYSCREEN);
  
  //FIXME : NOT tested !
  fWindow = ::CreateWindowEx(0, className,fName.c_str(),
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

  //G.Barrand : avoid to indirectly pass in
  //              WindowProc/[WM_SIZE,WM_PAINT]/This->DrawView()
  //            from this method. Else we have crash.
  fInCreateWindow = true;

  ::SetForegroundWindow(fWindow);
  ::ShowWindow(fWindow,SW_SHOWDEFAULT);
  ::UpdateWindow(fWindow);
  ::DrawMenuBar(fWindow);

  fInCreateWindow = false;
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLWin32Viewer::G4OpenGLWin32Viewer (
 G4OpenGLSceneHandler& scene
)
:G4VViewer (scene, -1)
,G4OpenGLViewer (scene)
,fMouseHovered(false)
,fMousePressed(false)
,fMousePressedX(0)
,fMousePressedY(0)
,fHDC(0)
,fWindow(0)
,fHGLRC(0)
,fInCreateWindow(false)
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
LRESULT CALLBACK G4OpenGLWin32Viewer::WindowProc(
 HWND aWindow
,UINT aMessage
,WPARAM aWParam
,LPARAM aLParam
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  switch (aMessage) {
    case WM_SIZE: {
      //FIXME : have to handle WM_RESIZE
      // Seems to be done (ovidio.pena AT upm.es, 2021/02/23)
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      if (This) {
        This->fWinSize_x = (G4int) LOWORD(aLParam);
        This->fWinSize_y = (G4int) HIWORD(aLParam);
        if (!This->fInCreateWindow) {
          This->SetView();
          glViewport(0, 0, This->fWinSize_x, This->fWinSize_y);
          This->DrawView();
        }
      }
      return 0;
    }

    case WM_PAINT: {
      PAINTSTRUCT ps;
      BeginPaint(aWindow, &ps);
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      if (This) {
        //FIXME : To have an automatic refresh someone have to redraw here.
        // Seems to be done (ovidio.pena AT upm.es, 2021/02/23)
        if(!This->fInCreateWindow) {
          This->SetView();
          This->ClearView();
          This->DrawView();
        }
      }
      EndPaint(aWindow, &ps);
      return 0;
    }

    case WM_LBUTTONDOWN: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      This->TrackMouse(LOWORD(aLParam), HIWORD(aLParam));
      return 0;
    }

    case WM_RBUTTONDOWN: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      This->TrackMouse(LOWORD(aLParam), HIWORD(aLParam));
      return 0;
    }

    case WM_LBUTTONUP: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      This->ReleaseMouse();
      return 0;
    }

    case WM_RBUTTONUP: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      This->ReleaseMouse();
      return 0;
    }

    case WM_MOUSEHOVER: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      This->fMouseHovered = true;
            return 0;
    }

    case WM_MOUSELEAVE: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);
      This->fMouseHovered = false;
      return 0;
    }

    case WM_MOUSEMOVE: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);

      if (!This->fMouseHovered) {
        // mouse hover/leave tracking
        TRACKMOUSEEVENT tme;
        tme.cbSize = sizeof(tme);
        tme.dwFlags = TME_HOVER | TME_LEAVE;
        tme.hwndTrack = aWindow;
        tme.dwHoverTime = HOVER_DEFAULT;
        ::TrackMouseEvent(&tme);
        This->fMouseHovered = true;
      }

      if (This->fMousePressed) {
        G4int x = (G4int) LOWORD(aLParam);
        G4int y = (G4int) HIWORD(aLParam);
        G4int dx = x - This->fMousePressedX;
        G4int dy = y - This->fMousePressedY;
        This->fMousePressedX = x;
        This->fMousePressedY = y;

        if (aWParam == MK_LBUTTON) {  // Rotation
          This->SetRotation(dx, dy);
        }

        if (aWParam == MK_RBUTTON) {  // Shift
          This->SetShift(dx, dy);
        }

        This->SetView();
        This->ClearView();
        This->DrawView();
      }

      return 0;
    }

    case WM_MOUSEWHEEL: {
      auto* This = (G4OpenGLWin32Viewer*)
                   ::GetWindowLongPtr(aWindow, GWLP_USERDATA);

      G4int delta = (short) HIWORD(aWParam);

      This->SetZoom(delta);

      This->SetView();
      This->ClearView();
      This->DrawView();
      return 0;
    }

    default:
      return DefWindowProc(aWindow, aMessage, aWParam, aLParam);
  }
//  return DefWindowProc(aWindow,aMessage,aWParam,aLParam);
}

//////////////////////////////////////////////////////////////////////////////
G4bool G4OpenGLWin32Viewer::SetWindowPixelFormat(
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
  
  G4int pixelIndex = ::ChoosePixelFormat(aHdc,&pfd);
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

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::TrackMouse(
 G4int x
,G4int y
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fMousePressed = true;
  fMousePressedX = x;
  fMousePressedY = y;
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::ReleaseMouse(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fMousePressed = false;
  fMousePressedX = 0;
  fMousePressedY = 0;
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::SetShift(
 G4int dx
,G4int dy
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  const G4double sceneRadius = GetSceneHandler()->GetScene()
                             ->GetExtent().GetExtentRadius();
  const G4double scale = 300;  // Roughly pixels per window, empirically chosen
  const G4double dxScene = dx*sceneRadius/scale;
  const G4double dyScene = dy*sceneRadius/scale;
  fVP.IncrementPan(-dxScene,dyScene);
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::SetRotation(
 G4int dx
,G4int dy
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  // Simple ad-hoc algorithms (borrowed from G4Qt3DViewer)
  const G4Vector3D& x_prime = fVP.GetViewpointDirection()
                              .cross(fVP.GetUpVector());
  const G4Vector3D& y_prime = x_prime.cross(fVP.GetViewpointDirection());
  const G4double scale = 200;  // Roughly pixels per window, empirically chosen
  G4Vector3D newViewpointDirection = fVP.GetViewpointDirection();
  newViewpointDirection += dx*x_prime/scale;
  newViewpointDirection += dy*y_prime/scale;
  fVP.SetViewpointDirection(newViewpointDirection.unit());

  if (fVP.GetRotationStyle() == G4ViewParameters::freeRotation) {
      G4Vector3D newUpVector = fVP.GetUpVector();
      newUpVector += dx*x_prime/scale;
      newUpVector += dy*y_prime/scale;
      fVP.SetUpVector(newUpVector.unit());
  }
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWin32Viewer::SetZoom(
 G4int delta
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
    if (fVP.GetFieldHalfAngle() == 0.) {  // Orthographic projection
        const G4double scale = 500;  // Empirically chosen
        fVP.MultiplyZoomFactor(1. + delta/scale);
    } else {                              // Perspective projection
        const G4double scale = fVP.GetFieldHalfAngle()/(10.*deg);  // Empirical
        fVP.SetDolly(fVP.GetDolly() + delta/scale);
    }
}

