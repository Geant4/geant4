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
// $Id: G4OpenInventorWinViewer.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
/*
 * jck : 05 Feb 1997 : Initial Implementation
 * jck : 21 Apr 1997 : Mods for SoXtHepViewer
 * gb : on Win32 use an SoXtExaminerViewer.
 * gb : 05 April 2004 : creation.
 * gb : 09 November 2004 : Pulldown menu with the escape menu item.
 * gb 14 November 2004 : inherit G4OpenInventorViewer.
 */

#ifdef G4VIS_BUILD_OIWIN32_DRIVER

// this :
#include "G4OpenInventorWinViewer.hh"

#include <Inventor/nodes/SoSelection.h>

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

#include "HEPVis/actions/SoGL2PSAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"
#include "G4VisManager.hh"

#include <windowsx.h>

// To have sizeChanged public :
class Geant4_SoWinExaminerViewer : public SoWinExaminerViewer {
public:
  Geant4_SoWinExaminerViewer(HWND parent,const char* name,SbBool embed)
  :SoWinExaminerViewer(parent,name,embed){}
  virtual void sizeChanged(const SbVec2s & size){
    SoWinExaminerViewer::sizeChanged(size);
  }
};

// File : 
#define ID_FILE_POSTSCRIPT 1
#define ID_FILE_PIXMAP_POSTSCRIPT 2
#define ID_FILE_INVENTOR 3
#define ID_FILE_ESCAPE 4
// Etc : 
#define ID_ETC_ERASE_DETECTOR 101
#define ID_ETC_ERASE_EVENT 102
#define ID_ETC_SET_SOLID 103
#define ID_ETC_SET_WIRE_FRAME 104
#define ID_ETC_SET_REDUCED_WIRE_FRAME 105
#define ID_ETC_SET_FULL_WIRE_FRAME 106
#define ID_ETC_SET_PREVIEW 107
#define ID_ETC_SET_PREVIEW_AND_FULL 108
#define ID_ETC_UPDATE_SCENE 109
#define ID_ETC_STATS 110
// Help :
#define ID_HELP_CONTROLS 201

//static void SecondaryLoopPostAction ();

static const char className[] = "G4OpenInventorShellWindow";

G4OpenInventorWinViewer::G4OpenInventorWinViewer(
 G4OpenInventorSceneHandler& sceneHandler
,const G4String& name)
:G4OpenInventorViewer (sceneHandler, name)
,fShell(0)
,fViewer(0)
{
  if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "Window name: " << fName << G4endl;
}


void G4OpenInventorWinViewer::Initialise() {

  G4String wName = fName;

  int width = 600;
  int height = 600;

  HWND parent = (HWND)fInteractorManager->GetParentInteractor ();
  if(!parent) {
    //Create a shell window :
    G4String shellName = wName;
    shellName += "_shell"; 
    static SbBool done = FALSE;
    if(done==FALSE) {
      HBRUSH brush = (HBRUSH) GetSysColorBrush(COLOR_BTNFACE);
      WNDCLASS wc;
      wc.style = CS_HREDRAW | CS_VREDRAW;
      wc.lpfnWndProc = (WNDPROC)WindowProc;
      wc.cbClsExtra = 0;
      wc.cbWndExtra = 0;
      wc.hInstance = ::GetModuleHandle(0);
      wc.hIcon = ::LoadIcon(0, IDI_APPLICATION);
      wc.hCursor = ::LoadCursor(0, IDC_ARROW);
      wc.hbrBackground = brush;
      wc.lpszMenuName = className;
      wc.lpszClassName = className;
      ::RegisterClass(&wc);
      done = TRUE;
    }
    
    width = fVP.GetWindowSizeHintX();
    height = fVP.GetWindowSizeHintX();
    
    HMENU menuBar = CreateMenu();

   {HMENU casc = CreatePopupMenu();
    ::AppendMenu(menuBar,MF_POPUP,(UINT)casc,"File");
    ::AppendMenu(casc,MF_STRING,ID_FILE_POSTSCRIPT,"PS (gl2ps)");
    ::AppendMenu(casc,MF_STRING,ID_FILE_PIXMAP_POSTSCRIPT,"PS (pixmap)");
    ::AppendMenu(casc,MF_STRING,ID_FILE_INVENTOR,"IV");
    ::AppendMenu(casc,MF_STRING,ID_FILE_ESCAPE,"Escape");}

   {HMENU casc = CreatePopupMenu();
    ::AppendMenu(menuBar,MF_POPUP,(UINT)casc,"Etc");
    ::AppendMenu(casc,MF_STRING,ID_ETC_ERASE_DETECTOR,"Erase detector");
    ::AppendMenu(casc,MF_STRING,ID_ETC_ERASE_EVENT,"Erase event");
    ::AppendMenu(casc,MF_STRING,ID_ETC_SET_SOLID,"Set solid");
    //::AppendMenu(casc,MF_STRING,ID_ETC_SET_WIRE_FRAME,"Set (G4) wire frame");
    ::AppendMenu(casc,MF_STRING,ID_ETC_SET_REDUCED_WIRE_FRAME,
                      "Set (G4) reduced wire frame");
    ::AppendMenu(casc,MF_STRING,ID_ETC_SET_FULL_WIRE_FRAME,
                      "Set (G4) full wire frame");
    ::AppendMenu(casc,MF_STRING,ID_ETC_SET_PREVIEW,
                      "Visible mothers + invisible daughters");
    ::AppendMenu(casc,MF_STRING,ID_ETC_SET_PREVIEW_AND_FULL,
                      "Visible mothers + visible daughters");
    ::AppendMenu(casc,MF_STRING,ID_ETC_UPDATE_SCENE,"Update scene");
    ::AppendMenu(casc,MF_STRING,ID_ETC_STATS,"Scene graph stats");}

   {HMENU casc = CreatePopupMenu();
    ::AppendMenu(menuBar,MF_POPUP,(UINT)casc,"Help");
    ::AppendMenu(casc,MF_STRING,ID_HELP_CONTROLS,"Controls");}

    fShell = ::CreateWindow(className, shellName.c_str(), 
                            WS_OVERLAPPEDWINDOW |
                            WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN,
                            CW_USEDEFAULT, CW_USEDEFAULT, 
                            width,height,
                            0,menuBar,::GetModuleHandle(0),0);
    // Retreive window and client sizez :
    RECT wrect,crect;
    GetWindowRect((HWND)fShell,&wrect);
    GetClientRect((HWND)fShell,&crect);
    int ww = wrect.right-wrect.left;
    int wh = wrect.bottom-wrect.top;
    int cw = crect.right-crect.left;
    int ch = crect.bottom-crect.top;
    // Compell client rect to be width height :
    MoveWindow((HWND)fShell,wrect.left,wrect.top,width+ww-cw,height+wh-ch,TRUE);
    ::SetWindowLongPtr((HWND)fShell,GWLP_USERDATA,LONG(this));
    ::SetWindowText((HWND)fShell,shellName.c_str());
    parent = fShell;
    fInteractorManager->AddShell(fShell);
  } else {
    char* str = fInteractorManager->GetCreationString();
    if(str!=0) wName = str;
  }
  fViewer = new Geant4_SoWinExaminerViewer(parent,wName.c_str(),TRUE);

  // Have a GL2PS render action :
  const SbViewportRegion& vpRegion = fViewer->getViewportRegion();
  fGL2PSAction = new SoGL2PSAction(vpRegion);
  fViewer->setGLRenderAction(fGL2PSAction);

  fViewer->setSize(SbVec2s(width,height));
  fViewer->setSceneGraph(fSoSelection);
  fViewer->viewAll();
  fViewer->saveHomePosition();
  fViewer->setTitle(fName);
  fViewer->show();
  if(fShell) {
    SoWin::show(fShell);
    fInteractorManager->FlushAndWaitExecution ();
  }
  fInteractorManager->SetCreatedInteractor (fViewer -> getWidget());
  fViewer->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_ADD);
}

G4OpenInventorWinViewer::~G4OpenInventorWinViewer () {
  if(fShell) fInteractorManager->RemoveShell(fShell);
  if(fViewer) {
    fViewer->setSceneGraph(0);
    delete fViewer;
  }
  if(fShell) {
    ::SetWindowLongPtr((HWND)fShell,GWLP_USERDATA,LONG(0));
    ::DestroyWindow((HWND)fShell);
  }
}

void G4OpenInventorWinViewer::FinishView () {
  if(!fViewer) return;
  fViewer->viewAll();
  fViewer->saveHomePosition();
}

void G4OpenInventorWinViewer::SetView () {
  G4OpenInventorViewer::SetView ();
  if(!fViewer) return;
  // Background.
  G4Colour b = fVP.GetBackgroundColour ();
  fViewer->setBackgroundColor
    (SbColor((float)b.GetRed(),(float)b.GetGreen(),(float)b.GetBlue()));
}
void G4OpenInventorWinViewer::ViewerRender () {
  if(!fViewer) return;
  fViewer->render();
}

SoCamera* G4OpenInventorWinViewer::GetCamera () {
  if(!fViewer) return 0;
  return fViewer->getCamera();
}


//////////////////////////////////////////////////////////////////////////////
LRESULT CALLBACK G4OpenInventorWinViewer::WindowProc ( 
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
    G4OpenInventorWinViewer* This = 
      (G4OpenInventorWinViewer*)::GetWindowLongPtr(aWindow,GWLP_USERDATA);
    if(This && This->fViewer) {
      This->fViewer->sizeChanged(SbVec2s(width,height));
    }
  }return 0;
  case WM_SETFOCUS:{ // Assume one child window !
    HWND hwnd = ::GetFirstChild(aWindow);
    if(hwnd!=0) ::SetFocus(hwnd);
  }return 0;
  case WM_DESTROY:{
    //G4OpenInventorWinViewer* This = 
    //  (G4OpenInventorWinViewer*)::GetWindowLongPtr(aWindow,GWLP_USERDATA);
    //::PostQuitMessage(0);
  }return 0;
  case WM_COMMAND:{
    G4OpenInventorWinViewer* This = 
      (G4OpenInventorWinViewer*)::GetWindowLongPtr(aWindow,GWLP_USERDATA);
    if(This) {
      if(aLParam==0) { //From menu.
        // File :
        if(aWParam==ID_FILE_POSTSCRIPT) {
          This->WritePostScript();
        } else if(aWParam==ID_FILE_PIXMAP_POSTSCRIPT) {
          This->WritePixmapPostScript();
        } else if(aWParam==ID_FILE_INVENTOR) {
          This->WriteInventor();
        } else if(aWParam==ID_FILE_ESCAPE) {
          This->Escape();
        // Etc :
        } else if(aWParam==ID_ETC_ERASE_DETECTOR) {
          This->EraseDetector();
        } else if(aWParam==ID_ETC_ERASE_EVENT) {
          This->EraseEvent();
        } else if(aWParam==ID_ETC_SET_SOLID) {
          This->SetSolid();
        } else if(aWParam==ID_ETC_SET_WIRE_FRAME) {
          This->SetWireFrame();
        } else if(aWParam==ID_ETC_SET_REDUCED_WIRE_FRAME) {
          This->SetReducedWireFrame(true);
        } else if(aWParam==ID_ETC_SET_FULL_WIRE_FRAME) {
          This->SetReducedWireFrame(false);
        } else if(aWParam==ID_ETC_SET_PREVIEW) {
          This->SetPreview();
        } else if(aWParam==ID_ETC_SET_PREVIEW_AND_FULL) {
          This->SetPreviewAndFull();
        } else if(aWParam==ID_ETC_UPDATE_SCENE) {
          This->UpdateScene();
        } else if(aWParam==ID_ETC_STATS) {
          This->SceneGraphStatistics();
        // Help :
        } else if(aWParam==ID_HELP_CONTROLS) {
          G4cout << This->Help() << G4endl;
        }
      }
    }
  }return 0;
  default:
    return (::DefWindowProc(aWindow,aMessage,aWParam,aLParam));
  }
}

#endif
