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
// $Id: G4OpenInventorXtViewer.cc,v 1.5 2004-11-09 09:16:47 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*
 * jck 05 Feb 1997 - Initial Implementation
 * jck 21 Apr 1997 
 *	Mods for SoXtHepViewer
 * gb : on Win32 use an SoXtExaminerViewer.
 * gb 05 April 2004 : revisit to separate Windows things.
 * gb 09 November 2004 : restore the escape button.
 * gb 09 November 2004 : have a menu bar in the viewer shell.
 * gb 09 November 2004 : have gl2ps file production.
 */
#ifdef G4VIS_BUILD_OIX_DRIVER

// this :
#include "G4OpenInventorXtViewer.hh"

#include <Inventor/nodes/SoSelection.h>

#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>

#include <X11/StringDefs.h>
#include <X11/Shell.h>

#include <Xm/Xm.h>
#include <Xm/PushB.h>
#include <Xm/Form.h>
#include <Xm/CascadeB.h>
#include <Xm/RowColumn.h>

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"

#include "HEPVis/actions/SoGL2PSAction.h"

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

  // Selection
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
    Arg args[10];
    char s[32];
    sprintf(s,"%dx%d",SIZE,SIZE);
    XtSetArg(args[0],XtNgeometry,XtNewString(s));
    XtSetArg(args[1],XtNborderWidth,0);
    XtSetArg(args[2],XtNtitle,XtNewString(wName.c_str()));
    fShell = XtAppCreateShell(shellName.c_str(),"Inventor",
	  		       topLevelShellWidgetClass,
			       SoXt::getDisplay(),
			       args,3); 

    XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
    XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
    XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
    XtSetArg(args[3],XmNbottomAttachment,XmATTACH_FORM);
    Widget form = XmCreateForm (fShell,(char*)"form",args,4);
    XtManageChild (form);

    XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
    XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
    XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
    Widget menuBar = XmCreateMenuBar (form,(char*)"menuBar",args,3);
    XtManageChild(menuBar);

    Widget menu = AddMenu(menuBar,"File","File");
    AddButton(menu,"PostScript",PostScriptButtonCbk);
    AddButton(menu,"Escape",EscapeButtonCbk);

    fViewer = new SoXtExaminerViewer(form,wName.c_str(),TRUE);
    
    XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_WIDGET);
    XtSetArg(args[1],XmNtopWidget       ,menuBar);
    XtSetArg(args[2],XmNleftAttachment  ,XmATTACH_FORM);
    XtSetArg(args[3],XmNrightAttachment ,XmATTACH_FORM);
    XtSetArg(args[4],XmNbottomAttachment,XmATTACH_FORM);
    XtSetValues(fViewer->getWidget(),args,5);

    fInteractorManager->AddShell(fShell);

  } else {
    char* str = fInteractorManager->GetCreationString();
    if(str!=0) wName = str;
    fViewer = new SoXtExaminerViewer(parent,wName.c_str(),TRUE);
  }

  // Have a GL2PS render action :
  const SbViewportRegion& vpRegion = fViewer->getViewportRegion();
  SoGL2PSAction* action = new SoGL2PSAction(vpRegion);
  fViewer->setGLRenderAction(action);

  // Else :
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
  if(fViewer) {
    fViewer->setSceneGraph(0);
    //FIXME : SGI : the below "delete" block things.
    //FIXME : CoinXt : the below "delete" crashe in ~SoXtRenderArea.
    //FIXME : delete fViewer;
  }
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

void G4OpenInventorXtViewer::WritePostScript(const G4String& aFile) {
  if(!fViewer) return;
  SoGL2PSAction* action = (SoGL2PSAction*)fViewer->getGLRenderAction();
  action->setFileName(aFile.c_str());
  action->enableFileWriting();
  fViewer->render();
  action->disableFileWriting();
  fViewer->render();
}

void G4OpenInventorXtViewer::EscapeButtonCbk(
 Widget,XtPointer aData,XtPointer) {
 G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
 This->fInteractorManager->RequireExitSecondaryLoop (OIV_EXIT_CODE);
}

void G4OpenInventorXtViewer::PostScriptButtonCbk(
 Widget,XtPointer aData,XtPointer) {
 G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
 This->WritePostScript();
}

Widget G4OpenInventorXtViewer::AddMenu(
 Widget aMenuBar
,const G4String& aName
,const G4String& aLabel
)
{
  // Pulldown menu :
  Widget menu = XmCreatePulldownMenu(aMenuBar,(char*)aName.c_str(),NULL,0);
  // Cascade button :
  Arg args[2];
  XmString cps = 
    XmStringLtoRCreate((char*)aLabel.c_str(),XmSTRING_DEFAULT_CHARSET);
  XtSetArg (args[0],XmNlabelString,cps);
  XtSetArg (args[1],XmNsubMenuId,menu);
  Widget widget = XmCreateCascadeButton(aMenuBar,(char*)aName.c_str(),args,2);
  XmStringFree (cps);
  XtManageChild(widget);
  return menu;
}
void G4OpenInventorXtViewer::AddButton (
 Widget aMenu
,const G4String& aLabel
,XtCallbackProc aCallback
)
{
  Widget widget = XmCreatePushButton(aMenu,(char*)aLabel.c_str(),NULL,0);
  XtManageChild(widget);
  XtAddCallback(widget,XmNactivateCallback,aCallback,(XtPointer)this);
}


#endif
