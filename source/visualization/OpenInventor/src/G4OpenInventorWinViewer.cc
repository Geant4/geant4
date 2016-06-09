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
// $Id: G4OpenInventorWinViewer.cc,v 1.2 2004/04/08 10:49:57 gbarrand Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
/*
 * jck 05 Feb 1997 - Initial Implementation
 * jck 21 Apr 1997 
 *	Mods for SoXtHepViewer
 * gb : on Win32 use an SoXtExaminerViewer.
 * gb 05 April 2004 : creation.
 */
#ifdef G4VIS_BUILD_OIWIN32_DRIVER

// this :
#include "G4OpenInventorWinViewer.hh"

#include <Inventor/nodes/SoSelection.h>

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"

#include <windowsx.h>

//
// Global variables 
//

//static void SecondaryLoopPostAction ();

static LRESULT CALLBACK WindowProc (HWND,UINT,WPARAM,LPARAM);

static const char className[] = "G4OpenInventorShellWindow";

void G4OpenInventorWinViewer::FinishView () {
  if(!fViewer) return;
  fViewer->viewAll();
  fViewer->saveHomePosition();
}

void G4OpenInventorWinViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  if (
      //??fG4OpenInventorSceneHandler.fPODLList.size() == 0 ||
      // We need a test for empty scene graph, such as
      // staticRoot.size() or something??????????  See temporary fix
      // in contructor.  (John Allison Aug 2001)
      CompareForKernelVisit(fG4OpenInventorSceneHandler.fLastVP)  ||
      CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();
  }      
  fLastVP = fVP;
  fG4OpenInventorSceneHandler.fLastVP = fVP;
}
 
G4bool G4OpenInventorWinViewer::CompareForKernelVisit
(G4ViewParameters&) {

  if (
      (fLastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (fLastVP.GetRepStyle ()        != fVP.GetRepStyle ())        ||
      (fLastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (fLastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (fLastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (fLastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (fLastVP.IsSection ()          != fVP.IsSection ())          ||
      // No need to visit kernel if section plane changes.
      (fLastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      (fLastVP.GetCutawayPlanes ().size () !=
                                 fVP.GetCutawayPlanes ().size ()) ||
      // No need to visit kernel if cutaway planes change.
      (fLastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (fLastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())
      ) {
      return true;;
  }
  if (fLastVP.IsDensityCulling () &&
      (fLastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  if (fLastVP.IsExplode () &&
      (fLastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;
      
  return false;
}

G4OpenInventorWinViewer::G4OpenInventorWinViewer
(G4OpenInventorSceneHandler& sceneHandler,
 const G4String& name)
:G4VViewer (sceneHandler, sceneHandler.IncrementViewCount(), name)
,fG4OpenInventorSceneHandler(sceneHandler)
,fShell(0)
,fViewer(0)
,fSelection(0)
,fInteractorManager(0)
{
  fNeedKernelVisit = true;  //?? Temporary, until KernelVisitDecision fixed.

  fInteractorManager = 
    ((G4OpenInventor*)fG4OpenInventorSceneHandler.GetGraphicsSystem())->
    GetInteractorManager();
  //Widget toplevel = (Widget)fInteractorManager->GetMainInteractor ();

  //fInteractorManager->
  //AddSecondaryLoopPostAction((G4SecondaryLoopAction)SecondaryLoopPostAction);

  G4cout << "Window name: " << fName << G4endl;
  // 
  // Selection
  //
  fSelection = new SoSelection;
  fSelection->policy = SoSelection::SINGLE;
  fSelection->ref();
  fSelection->addChild(fG4OpenInventorSceneHandler.root);

  G4String wName = fName;

#define SIZE 400
  HWND parent = (HWND)fInteractorManager->GetParentInteractor ();
  if(!parent) {
    //Create a shell window :
    G4String shellName = wName;
    shellName += "_shell"; 
    static SbBool done = FALSE;
    if(done==FALSE) {
      WNDCLASS wc;
      wc.style = CS_HREDRAW | CS_VREDRAW;
      wc.lpfnWndProc = (WNDPROC)WindowProc;
      wc.cbClsExtra = 0;
      wc.cbWndExtra = 0;
      wc.hInstance = ::GetModuleHandle(0);
      wc.hIcon = ::LoadIcon(0, IDI_APPLICATION);
      wc.hCursor = ::LoadCursor(0, IDC_ARROW);
      wc.hbrBackground = 0;
      wc.lpszMenuName = className;
      wc.lpszClassName = className;
      ::RegisterClass(&wc);
      done = TRUE;
    }
    //  Compell window to be created at 0,0 to bypass 
    // the 'black border' problem. 
    fShell = ::CreateWindow(className, shellName.c_str(), 
                            WS_OVERLAPPEDWINDOW,
                            //CW_USEDEFAULT, CW_USEDEFAULT, 
                            0,0,400,400,0, 0,::GetModuleHandle(0),0);
    // Retreive window and client sizez :
    RECT wrect,crect;
    GetWindowRect((HWND)fShell,&wrect);
    GetClientRect((HWND)fShell,&crect);
    int ww = wrect.right-wrect.left;
    int wh = wrect.bottom-wrect.top;
    int cw = crect.right-crect.left;
    int ch = crect.bottom-crect.top;
    // Compell client rect to be 400 400 :
    MoveWindow((HWND)fShell,wrect.left,wrect.top,400+ww-cw,400+wh-ch,TRUE);
    ::SetWindowLong((HWND)fShell,GWL_USERDATA,LONG(this));
    ::SetWindowText((HWND)fShell,shellName.c_str());
    parent = fShell;
    fInteractorManager->AddShell(fShell);
  } else {
    char* str = fInteractorManager->GetCreationString();
    if(str!=0) wName = str;
  }
  fViewer = new SoWinExaminerViewer(parent,wName.c_str(),TRUE);
  fViewer->setSize(SbVec2s(SIZE,SIZE));
  fViewer->setSceneGraph(fSelection);
  fViewer->viewAll();
  fViewer->saveHomePosition();
  fViewer->setTitle(fName);
  fViewer->show();
  if(fShell) {
    SoWin::show(fShell);
    fInteractorManager->FlushAndWaitExecution ();
  }
  fInteractorManager->SetCreatedInteractor (fViewer -> getWidget());
}

G4OpenInventorWinViewer::~G4OpenInventorWinViewer () {
  if(fShell) fInteractorManager->RemoveShell(fShell);
  if(fViewer) delete fViewer;
  if(fShell) {
    ::SetWindowLong((HWND)fShell,GWL_USERDATA,LONG(0));
    ::DestroyWindow((HWND)fShell);
  }
  if(fSelection) fSelection->unref();
}

void G4OpenInventorWinViewer::ClearView () {
}

void G4OpenInventorWinViewer::SetView () {
}

void G4OpenInventorWinViewer::DrawView () {
  G4cout << "debug Iv::DrawViewer " <<G4endl;
  KernelVisitDecision();
  ProcessView();
  FinishView();
}

void G4OpenInventorWinViewer::ShowView () {
  fInteractorManager -> SecondaryLoop ();
}
//////////////////////////////////////////////////////////////////////////////
LRESULT CALLBACK WindowProc ( 
 HWND   aWindow
,UINT   aMessage
,WPARAM aWParam
,LPARAM aLParam
)
//////////////////////////////////////////////////////////////////////////////
//  Below treatment of WM_SIZE, WM_SETFOCUS not necessary 
// with TGS, but needed with SoFree. WM_DESTROY needed for 
// 'main top level window' so that 'Close window' induces
// the end of the task.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  switch (aMessage) { 
  case WM_SIZE:{ // Assume one child window !
    int width = LOWORD(aLParam);
    int height = HIWORD(aLParam);
    //printf("debug : G4SoWindow : WMS_SIZE : %d %d\n",width,height);
    HWND hwnd = ::GetFirstChild(aWindow);
    if(hwnd!=0) {
      ::MoveWindow(hwnd,0,0,width,height,TRUE);
    }
  }return 0;
  case WM_SETFOCUS:{ // Assume one child window !
    HWND hwnd = ::GetFirstChild(aWindow);
    if(hwnd!=0) ::SetFocus(hwnd);
  }return 0;
  case WM_DESTROY:{
    //G4OpenInventorWinViewer* This = 
    //  (G4OpenInventorWinViewer*)::GetWindowLong(aWindow,GWL_USERDATA);
    //::PostQuitMessage(0);
  }return 0;
  default:
    return (::DefWindowProc(aWindow,aMessage,aWParam,aLParam));
  }
}

#endif
