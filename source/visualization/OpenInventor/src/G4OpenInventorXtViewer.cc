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
// $Id: G4OpenInventorXtViewer.cc,v 1.13 2004-11-15 14:53:30 gbarrand Exp $
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
 * gb 14 November 2004 : inherit G4OpenInventorViewer.
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

#include "HEPVis/actions/SoGL2PSAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"

G4OpenInventorXtViewer::G4OpenInventorXtViewer(
 G4OpenInventorSceneHandler& sceneHandler
,const G4String& name)
:G4OpenInventorViewer (sceneHandler, name)
,fShell(0)
,fViewer(0)
{
  G4cout << "Window name: " << fName << G4endl;

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

   {Widget menu = AddMenu(menuBar,"File","File");
    AddButton(menu,"PS (gl2ps)",PostScriptCbk);
    AddButton(menu,"PS (pixmap)",PixmapPostScriptCbk);
    AddButton(menu,"IV",WriteInventorCbk);
    AddButton(menu,"Escape",EscapeCbk);}

   {Widget menu = AddMenu(menuBar,"Etc","Etc");
    AddButton(menu,"Triangles",CountTrianglesCbk);
    AddButton(menu,"Erase detector",EraseDetectorCbk);
    AddButton(menu,"Erase event",EraseEventCbk);}

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

  fViewer->setSize(SbVec2s(SIZE,SIZE));

  // Have a GL2PS render action :
  const SbViewportRegion& vpRegion = fViewer->getViewportRegion();
  fGL2PSAction = new SoGL2PSAction(vpRegion);
  fViewer->setGLRenderAction(fGL2PSAction);

  // Else :
  fViewer->setSceneGraph(fSoSelection);
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
}

void G4OpenInventorXtViewer::FinishView () {
  if(!fViewer) return;
  fViewer->viewAll();
  fViewer->saveHomePosition();
}

void G4OpenInventorXtViewer::ViewerRender () {
  if(!fViewer) return;
  fViewer->render();
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

void G4OpenInventorXtViewer::EscapeCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->Escape();
}

void G4OpenInventorXtViewer::PostScriptCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->WritePostScript();
}

void G4OpenInventorXtViewer::PixmapPostScriptCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->WritePixmapPostScript();
}

void G4OpenInventorXtViewer::CountTrianglesCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->CountTriangles();
}

void G4OpenInventorXtViewer::WriteInventorCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->WriteInventor();
}

void G4OpenInventorXtViewer::EraseDetectorCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->EraseDetector();
}

void G4OpenInventorXtViewer::EraseEventCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtViewer* This = (G4OpenInventorXtViewer*)aData;
  This->EraseEvent();
}


#endif


