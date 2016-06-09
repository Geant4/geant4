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
// $Id: G4OpenInventorXtViewer.cc,v 1.1 2004/04/08 09:41:11 gbarrand Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
/*
 * jck 05 Feb 1997 - Initial Implementation
 * jck 21 Apr 1997 
 *	Mods for SoXtHepViewer
 * gb : on Win32 use an SoXtExaminerViewer.
 * gb 05 April 2004 : revisit to separate Windows things.
 */
#ifdef G4VIS_BUILD_OIX_DRIVER

// this :
#include "G4OpenInventorXtViewer.hh"

#include <Inventor/nodes/SoSelection.h>

#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <X11/StringDefs.h>

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"

//
// Global variables 
//

//static void SecondaryLoopPostAction ();

void G4OpenInventorXtViewer::FinishView () {
  if(!fViewer) return;
  fViewer->viewAll();
  fViewer->saveHomePosition();
}

//static void SecondaryLoopPostAction ()
//{
//  Display *display = fViewer->getDisplay();
//  XSync(display, False);
//  if(interactorManager -> GetExitSecondaryLoopCode ()==OIV_EXIT_CODE) {
//    if (fShell!=0) XtRealizeWidget(fShell);
//  }
//}

void G4OpenInventorXtViewer::KernelVisitDecision () {
  
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
 
G4bool G4OpenInventorXtViewer::CompareForKernelVisit
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

G4OpenInventorXtViewer::G4OpenInventorXtViewer
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
  Widget parent = (Widget)fInteractorManager->GetParentInteractor ();
  if(!parent) {  
    //Create a shell window :
    G4String shellName = wName;
    shellName += "_shell"; 
    Arg args[3];
    char s[32];
    sprintf(s,"%dx%d",SIZE,SIZE);
    XtSetArg(args[0],XtNgeometry,XtNewString(s));
    XtSetArg(args[1],XtNborderWidth,0);
    XtSetArg(args[2],XtNtitle,XtNewString(wName.c_str()));
    fShell = XtAppCreateShell(shellName.c_str(),"Inventor",
	  		       topLevelShellWidgetClass,
			       SoXt::getDisplay(),
			       args,3); 
    parent = fShell;
    fInteractorManager->AddShell(fShell);
  } else {
    char* str = fInteractorManager->GetCreationString();
    if(str!=0) wName = str;
  }
  //
  // Create and Customize the Viewer
  //
  fViewer = new SoXtExaminerViewer(parent,wName.c_str(),TRUE);
  fViewer->setSize(SbVec2s(SIZE,SIZE));
  fViewer->setSceneGraph(fSelection);
  fViewer->viewAll();
  fViewer->saveHomePosition();
  fViewer->setTitle(fName);
  fViewer->show();
  if(fShell) {
    SoXt::show(fShell);
    fInteractorManager->FlushAndWaitExecution ();
  }
  fInteractorManager->SetCreatedInteractor (fViewer -> getWidget());
}

G4OpenInventorXtViewer::~G4OpenInventorXtViewer () {
  if(fShell) fInteractorManager->RemoveShell(fShell);
  if(fViewer) delete fViewer;
  if(fShell) XtDestroyWidget(fShell);
  if(fSelection) fSelection->unref();
}

void G4OpenInventorXtViewer::ClearView () {
}

void G4OpenInventorXtViewer::SetView () {
}

void G4OpenInventorXtViewer::DrawView () {
  G4cout << "debug Iv::DrawViewer " <<G4endl;
  KernelVisitDecision();
  ProcessView();
  FinishView();
}

void G4OpenInventorXtViewer::ShowView () {
  fInteractorManager -> SecondaryLoop ();
}

#endif
