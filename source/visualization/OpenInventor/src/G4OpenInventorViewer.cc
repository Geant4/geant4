// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorViewer.cc,v 1.4 2001-02-23 15:43:13 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*
 * jck 05 Feb 1997 - Initial Implementation
 * jck 21 Apr 1997 
 *	Mods for SoXtHepViewer
 * gb : on Win32 use an SoXtExaminerViewer.
 */
#ifdef G4VIS_BUILD_OI_DRIVER

#include <Inventor/Xt/SoXt.h>
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoSelection.h>

#include <HEPVis/viewers/SoWindow.h>

#ifdef WIN32
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#else
#include <HEPVis/Xt/viewers/SoXtHepViewer.h>
#endif

#include "G4OpenInventor.hh"
#include "G4OpenInventorViewer.hh"
#include "G4OpenInventorSceneHandler.hh"

#include "G4VInteractorManager.hh"

//
// Global variables 
//

//static void SecondaryLoopPostAction ();

void G4OpenInventorViewer::FinishView () {
  if(fViewer==NULL) return;
#ifdef WIN32
  ((SoXtViewer*)fViewer)->viewAll();
  ((SoXtViewer*)fViewer)->saveHomePosition();
#else
  ((SoXtHepViewer*)fViewer)->getCurrentViewer()->viewAll();
  ((SoXtHepViewer*)fViewer)->getCurrentViewer()->saveHomePosition();
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
//  Display *display = fViewer->getDisplay();
//  XSync(display, False);
//  if(interactorManager -> GetExitSecondaryLoopCode ()==OIV_EXIT_CODE) {
//    if (fShell!=NULL) XtRealizeWidget(fShell);
//  }
//}

void G4OpenInventorViewer::KernelVisitDecision () {
  
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
      (lastVP.GetCutawayPlanes ().size () !=
                                 fVP.GetCutawayPlanes ().size ()) ||

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

G4OpenInventorViewer::G4OpenInventorViewer (G4OpenInventorSceneHandler& scene,
					const G4String& name)
:G4VViewer (scene, scene.IncrementViewCount(), name)
,fSceneHandler(scene)
,fShell(NULL)
,fWindow(NULL)
,fViewer(NULL)
,fSelection(NULL)
,fInteractorManager(NULL)
{
  fInteractorManager = 
    ((G4OpenInventor*)fSceneHandler.GetGraphicsSystem())->
    GetInteractorManager();
  Widget toplevel = (Widget)fInteractorManager->GetMainInteractor ();

//fInteractorManager->
//AddSecondaryLoopPostAction ((G4SecondaryLoopAction)SecondaryLoopPostAction);

  G4cout << "Window name: " << fName << G4endl;
  // 
  // Selection
  //
  fSelection = new SoSelection;
  fSelection->policy = SoSelection::SINGLE;
  fSelection->ref();
  fSelection->addChild(fSceneHandler.root);

  Widget    parent = (Widget)fInteractorManager->GetParentInteractor ();

  G4String wName = fName;
  const char* wname = wName.data();

#define SIZE 400
  if(parent==NULL) {  //Create a shell window :
    fWindow = new SoWindow(wname);
    fWindow->setTitle(wname);
    fWindow->setSize(SbVec2s(SIZE,SIZE));
    fShell = parent = fWindow->getWidget();
    fInteractorManager->AddShell(fShell);
  } else {
    char* str = fInteractorManager->GetCreationString ();
    if(str!=NULL) wname = str;
  }
  //
  // Create and Customize the Viewer
  //
#ifdef WIN32
  fViewer = new SoXtExaminerViewer(parent,wname,TRUE);
  fViewer->setSize(SbVec2s(SIZE,SIZE));
  ((SoXtViewer*)fViewer)->setSceneGraph(fSelection);
  ((SoXtViewer*)fViewer)->viewAll();
  ((SoXtViewer*)fViewer)->saveHomePosition();
#else
  fViewer = new SoXtHepViewer (parent,wname,TRUE);
  ((SoXtHepViewer*)fViewer)->hideTextArea();
  ((SoXtHepViewer*)fViewer)->resizeMain(SIZE,SIZE);
  ((SoXtHepViewer*)fViewer)->setQuitCallback(quitCB, fInteractorManager);
  ((SoXtHepViewer*)fViewer)->setSceneGraph(fSelection);
  ((SoXtHepViewer*)fViewer)->getCurrentViewer()->viewAll();
  ((SoXtHepViewer*)fViewer)->getCurrentViewer()->saveHomePosition();
#endif
  fViewer->setTitle(fName);
  fViewer->show();
 
  if(fWindow!=NULL) {
    fWindow->show();
    fInteractorManager->FlushAndWaitExecution ();
  }

  fInteractorManager->SetCreatedInteractor (fViewer -> getWidget());
}

G4OpenInventorViewer::~G4OpenInventorViewer () {
  if(fShell!=NULL) fInteractorManager -> RemoveShell (fShell);
  if(fViewer!=NULL) {
#ifdef WIN32
    delete ((SoXtExaminerViewer*)fViewer);
#else
    delete ((SoXtHepViewer*)fViewer);
#endif
  }
  if(fWindow!=NULL) delete fWindow;
  if(fSelection!=NULL) fSelection->unref();
}

void G4OpenInventorViewer::ClearView () {
}

void G4OpenInventorViewer::SetView () {
}

void G4OpenInventorViewer::DrawView () {
  G4cout << "debug Iv::DrawViewer " <<G4endl;
  KernelVisitDecision ();
  ProcessView         ();
  FinishView          ();
}

void G4OpenInventorViewer::ShowView () {
  fInteractorManager -> SecondaryLoop ();
}

#endif








