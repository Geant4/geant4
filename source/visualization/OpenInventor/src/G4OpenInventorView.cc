// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorView.cc,v 1.1 1999-01-07 16:15:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*
 * jck 05 Feb 1997 - Initial Implementation
 * jck 21 Apr 1997 
 *	Mods for SoXtHepViewer
 *
 */
#ifdef G4VIS_BUILD_OI_DRIVER

#ifndef WIN32  //SoXtHepViewer not ported on Windows/NT.
#define HEP_VIEWER
#endif

#include <Inventor/Xt/SoXt.h>
#include <Inventor/nodes/SoSelection.h>

#ifdef HEP_VIEWER
#include <HEPVis/Xt/viewers/SoXtHepViewer.h>
#else
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#endif

#include "G4OpenInventor.hh"
#include "G4OpenInventorView.hh"
#include "G4OpenInventorScene.hh"

#include "G4VInteractorManager.hh"

//
// Global variables 
//

//static void SecondaryLoopPostAction ();

void G4OpenInventorView::FinishView () {
#ifdef HEP_VIEWER
  ((SoXtHepViewer*)G4OIViewer)->getCurrentViewer()->viewAll();
  ((SoXtHepViewer*)G4OIViewer)->getCurrentViewer()->saveHomePosition();
#else
  ((SoXtViewer*)G4OIViewer)->viewAll();
  ((SoXtViewer*)G4OIViewer)->saveHomePosition();
#endif
}

//
// Exit
//
int quitCB(void* interactorManager) {
  ((G4VInteractorManager*)interactorManager) -> 
    RequireExitSecondaryLoop (OIV_EXIT_CODE);
  return 0;
}

//static void SecondaryLoopPostAction ()
//{
//  Display *display = G4OIViewer->getDisplay();
//  XSync(display, False);
//  if(interactorManager -> GetExitSecondaryLoopCode ()==OIV_EXIT_CODE) {
//    if (G4OIShell!=NULL) XtRealizeWidget(G4OIShell);
//  }
//}

void G4OpenInventorView::KernelVisitDecision () {
  
  //
  // Trigger a display List refresh if necessary.  This is a checklist
  // of relevant view parameters.
  //
  static G4ViewParameters lastVP;  // Initialised to default.
  G4bool need = false;
  if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.GetRepStyle ()        != fVP.GetRepStyle ())        ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (lastVP.IsSection ()          != fVP.IsSection ())          ||

      (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      (lastVP.GetCutawayPlanes ().entries () !=
                              fVP.GetCutawayPlanes ().entries ()) ||

      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())
      ) {
      need = true;;
  }
  if (!need && lastVP.IsDensityCulling () &&
      (lastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    need = true;

  if (!need && lastVP.IsExplode () &&
      (lastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    need = true;
      
  if (need) {
    lastVP = fVP;
    NeedKernelVisit ();
  }
}

G4OpenInventorView::G4OpenInventorView (G4OpenInventorScene& scene,
					const G4String& name):
G4VView (scene, scene.IncrementViewCount (), name),
fScene (scene),
OIvisualfound (false),
G4OIViewer(NULL),
G4OIShell(NULL)
{
  interactorManager = ((G4OpenInventor*)fScene.GetGraphicsSystem()) -> 
    GetInteractorManager ();
  Widget toplevel = (Widget)interactorManager->GetMainInteractor ();

//interactorManager->
//AddSecondaryLoopPostAction ((G4SecondaryLoopAction)SecondaryLoopPostAction);

  OIvisualfound = true;

  G4cout << "Window name: " << fName << endl;
  // 
  // Selection
  //
  G4OISelection = new SoSelection;
  G4OISelection->policy = SoSelection::SINGLE;
  G4OISelection->ref();
  G4OISelection->addChild(fScene.root);

  Widget    parent = (Widget)interactorManager->GetParentInteractor ();
  SbBool    buildInsideParent;
  char charOIName[100];
  strncpy (charOIName, fName, 100);
  char*     wname = charOIName;

  if(parent==NULL) {  //Ask Inventor to create a viewer in its own shell.
    parent            = toplevel;
    buildInsideParent = FALSE;
  } else {            //Ask Inventor to create a viewer in the given parent.
    buildInsideParent = TRUE;
    char* str = interactorManager->GetCreationString ();
    if(str!=NULL) wname = str;
  }
  //
  // Create and Customize the Viewer
  //
#ifdef HEP_VIEWER
  G4OIViewer = new SoXtHepViewer (parent,wname,buildInsideParent);
  ((SoXtHepViewer*)G4OIViewer)->hideTextArea();
  ((SoXtHepViewer*)G4OIViewer)->resizeMain(550,550);
  ((SoXtHepViewer*)G4OIViewer)->setQuitCallback(quitCB, interactorManager);
  ((SoXtHepViewer*)G4OIViewer)->setSceneGraph(G4OISelection);
  ((SoXtHepViewer*)G4OIViewer)->getCurrentViewer()->viewAll();
  ((SoXtHepViewer*)G4OIViewer)->getCurrentViewer()->saveHomePosition();
#else
  G4OIViewer = new SoXtExaminerViewer(parent,wname,buildInsideParent);
  G4OIViewer->setSize(SbVec2s(550,550));
  ((SoXtViewer*)G4OIViewer)->setSceneGraph(G4OISelection);
  ((SoXtViewer*)G4OIViewer)->viewAll();
  ((SoXtViewer*)G4OIViewer)->saveHomePosition();
#endif
  G4OIViewer->setTitle(fName);
  G4OIViewer->show();
 
#if defined(WIN32) && !defined(SoFreePackage)
  //  buildInsideParent==FALSE, isTopLevel, getShellWidget 
  // do not work with Inventor-2.4 on Win32 !
  // All viewers are going to be put in the same shell window !
  if(buildInsideParent==FALSE) {
    G4OIShell  = toplevel; 
    interactorManager->AddShell (G4OIShell);
    SoXt::show(G4OIShell);
    interactorManager->FlushAndWaitExecution ();
  }
#else
  if(G4OIViewer->isTopLevelShell()==TRUE) {
    G4OIShell  = G4OIViewer -> getShellWidget();
    if(G4OIShell!=NULL) {
      interactorManager->AddShell (G4OIShell);
      SoXt::show(G4OIShell);
      interactorManager->FlushAndWaitExecution ();
    }
  }
#endif

  interactorManager->SetCreatedInteractor (G4OIViewer -> getWidget());
}

G4OpenInventorView::~G4OpenInventorView () {
  if(G4OIShell!=NULL) {
    interactorManager -> RemoveShell (G4OIShell);
  }
#ifdef HEP_VIEWER
  delete ((SoXtHepViewer*)G4OIViewer);
#else
  delete ((SoXtExaminerViewer*)G4OIViewer);
#endif
  if(G4OISelection!=NULL) {
    G4OISelection->unref();
  }
}

void G4OpenInventorView::ClearView () {
}

void G4OpenInventorView::SetView () {
}

void G4OpenInventorView::DrawView () {
  KernelVisitDecision ();
  ProcessView         ();
  FinishView          ();
}

void G4OpenInventorView::ShowView () {
  interactorManager -> SecondaryLoop ();
}

#endif








