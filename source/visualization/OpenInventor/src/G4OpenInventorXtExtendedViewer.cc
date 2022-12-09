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

// this :
#include "G4OpenInventorXtExtendedViewer.hh"

#include <Inventor/nodes/SoSelection.h>

#include <Inventor/Xt/SoXt.h>

//Replaces inclusion of SoXtExaminerViewer.h
#include <Inventor/Xt/viewers/SoXtFlyViewer.h>

#include <X11/StringDefs.h>
#include <X11/Shell.h>

#include <Xm/Xm.h>
#include <Xm/PushB.h>
#include <Xm/Form.h>
#include <Xm/CascadeB.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>

#include <cstdio>

#include "HEPVis/actions/SoGL2PSAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorXtExaminerViewer.hh"
#include "G4VInteractorManager.hh"
#include "G4VisManager.hh"
#include "G4AttCheck.hh"

G4OpenInventorXtExtendedViewer::G4OpenInventorXtExtendedViewer(
 G4OpenInventorSceneHandler& sceneHandler
,const G4String& name)
:G4OpenInventorViewer (sceneHandler, name)
,fShell(0)
,fViewer(0)
,fHelpForm(0)
,fHelpText(0)
{
  if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "Window name: " << fName << G4endl;
}


void G4OpenInventorXtExtendedViewer::Initialise() {
  
  G4String wName = fName;
  
  Widget parent = (Widget)fInteractorManager->GetParentInteractor ();
  // G4cout << "DEBUG G4OpenInventorXtExtendedViewer: parent = " << parent << G4endl;
  int width = 600;
  int height = 600;

  if(!parent) {  
    // Check if user has specified an X-Windows-type geometry string...
    char s[32];

    G4String sgeometry = fVP.GetXGeometryString();
    if(sgeometry.empty()) {
      G4cout << "ERROR: Geometry string \""
             << sgeometry	   
             << "\" is empty.  Using \"600x600\"."
             << G4endl;
      width = 600;
      height = 600;  
      snprintf(s,32,"%dx%d",width,height);
      sgeometry = s;
    } else {
      width = fVP.GetWindowSizeHintX();
      height = fVP.GetWindowSizeHintX();
    }

    //Create a shell window :
    G4String shellName = wName;
    shellName += "_shell"; 
    Arg args[10];
    XtSetArg(args[0],XtNgeometry,XtNewString(sgeometry.c_str()));
    XtSetArg(args[1],XtNborderWidth,0);
    XtSetArg(args[2],XtNtitle,XtNewString(wName.c_str()));
    fShell = XtAppCreateShell(shellName.c_str(),"Inventor",
	  		       topLevelShellWidgetClass,
			       SoXt::getDisplay(),
			       args,3);
    
    // G4cout << "DEBUG CREATING THE VIEWER WITH CREATED SHELL = " << fShell << G4endl;
    fViewer = new G4OpenInventorXtExaminerViewer(fShell, wName.c_str(), TRUE);
    fViewer->addEscapeCallback(EscapeFromKeyboardCbk, (void *)this);
     
    // FWJ (viewpoints don't work with this!)
    //    fViewer->setAutoClipping((SbBool)0);

    //XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
    //XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
    //XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
    //XtSetArg(args[3],XmNbottomAttachment,XmATTACH_FORM);
    //Widget form = XmCreateForm (fShell,(char*)"form",args,4);
    //XtManageChild (form);

    Widget menuBar = fViewer->getMenuBar();
    
    //XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
    //XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
    //XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
    //Widget menuBar = XmCreateMenuBar (form,(char*)"menuBar",args,3);
    //XtManageChild(menuBar);

   {Widget menu = fViewer->getMenu();
   //{Widget menu = AddMenu(menuBar,"File","File");
    AddButton(menu,"Write PS (gl2ps)",PostScriptCbk);
    AddButton(menu, "Write PDF (gl2ps)", PDFCbk);
    AddButton(menu,"Write PS (pixmap)",PixmapPostScriptCbk);
    AddButton(menu,"Write IV",WriteInventorCbk);
    AddButton(menu,"Escape",EscapeCbk);}

   {Widget menu = AddMenu(menuBar,"Etc","Etc");
    AddButton(menu,"Erase detector",EraseDetectorCbk);
    AddButton(menu,"Erase event",EraseEventCbk);
    AddButton(menu,"Set solid",SetSolidCbk);
/*    AddButton(menu,"Set (G4) wire frame",SetWireFrameCbk);*/
    AddButton(menu,"Set (G4) reduced wire frame",SetReducedWireFrameCbk);
    AddButton(menu,"Set (G4) full wire frame",SetFullWireFrameCbk);
    AddButton(menu,"Visible mothers + invisible daughters",SetPreviewCbk);
    AddButton(menu,"Visible mothers + visible daughters",SetPreviewAndFullCbk);
    AddButton(menu,"Update scene",UpdateSceneCbk);
    AddButton(menu,"Scene graph stats",SceneGraphStatisticsCbk);
   }

   {Widget menu = AddMenu(menuBar,"Help","Help");
    AddButton(menu,"Controls",HelpCbk);}

    //fViewer = new SoXtExaminerViewer(form,wName.c_str(),TRUE);
    
    XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_WIDGET);
    XtSetArg(args[1],XmNtopWidget       ,menuBar);
    XtSetArg(args[2],XmNleftAttachment  ,XmATTACH_FORM);
    XtSetArg(args[3],XmNrightAttachment ,XmATTACH_FORM);
    XtSetArg(args[4],XmNbottomAttachment,XmATTACH_FORM);
    XtSetValues(fViewer->getWidget(),args,5);

    fHelpForm = XmCreateFormDialog(fShell,(char*)"help",NULL,0);
    XtSetArg(args[0],XmNleftAttachment  ,XmATTACH_FORM);
    XtSetArg(args[1],XmNrightAttachment ,XmATTACH_FORM);
    XtSetArg(args[2],XmNbottomAttachment,XmATTACH_FORM);
    Widget cancel = XmCreatePushButton(fHelpForm,(char*)"helpCancel",args,3);
    XtAddCallback(cancel,XmNactivateCallback,HelpCancelCbk,(XtPointer)this);
    XtManageChild(cancel);
    XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
    XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
    XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
    XtSetArg(args[3],XmNbottomAttachment,XmATTACH_WIDGET);
    XtSetArg(args[4],XmNbottomWidget    ,cancel);
    fHelpText = XmCreateScrolledText(fHelpForm,(char*)"helpText",args,5);
    XtManageChild(fHelpText);

    fInteractorManager->AddShell(fShell);

  } else {
    char* str = fInteractorManager->GetCreationString();
    if(str!=0) wName = str;
    // G4cout << "DEBUG CREATING THE VIEWER WITH parent = " << parent << G4endl;
    fViewer = new G4OpenInventorXtExaminerViewer(parent, wName.c_str(), TRUE);
  }

  // Use our own SelectionCB for the Xt viewer to allow for abbreviated output
  // when picking a trajectory
  fSoSelection->removeSelectionCallback(G4OpenInventorViewer::SelectionCB,
                                        this);
//  fSoSelection->addSelectionCallback(SelectionCB, this);

  fViewer->setSize(SbVec2s(width,height));

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
  // TJR added:
  fViewer->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_ADD);
}

G4OpenInventorXtExtendedViewer::~G4OpenInventorXtExtendedViewer () {
  if(fShell) fInteractorManager->RemoveShell(fShell);
  if(fViewer) {
    fViewer->setSceneGraph(0);
    //FIXME : SGI : the below "delete" block things.
    //FIXME : CoinXt : the below "delete" crashe in ~SoXtRenderArea.
    //FIXME : delete fViewer;
  }
  if(fShell) XtDestroyWidget(fShell);
}

void G4OpenInventorXtExtendedViewer::FinishView () {
  if(!fViewer) return;
  fViewer->viewAll();
  fViewer->saveHomePosition();
}

void G4OpenInventorXtExtendedViewer::SetView () {
  G4OpenInventorViewer::SetView ();
  if(!fViewer) return;
  // Background.
  G4Colour b = fVP.GetBackgroundColour ();
  fViewer->setBackgroundColor
    (SbColor((float)b.GetRed(),(float)b.GetGreen(),(float)b.GetBlue()));
}


void G4OpenInventorXtExtendedViewer::ViewerRender () {
  if(!fViewer) return;
  fViewer->render();
}

SoCamera* G4OpenInventorXtExtendedViewer::GetCamera () {
  if(!fViewer) return 0;
  return fViewer->getCamera();
}

Widget G4OpenInventorXtExtendedViewer::AddMenu(
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
    XmStringLtoRCreate((char*)aLabel.c_str(),(char*)XmSTRING_DEFAULT_CHARSET);
  XtSetArg (args[0],XmNlabelString,cps);
  XtSetArg (args[1],XmNsubMenuId,menu);
  Widget widget = XmCreateCascadeButton(aMenuBar,(char*)aName.c_str(),args,2);
  XmStringFree (cps);
  XtManageChild(widget);
  return menu;
}
void G4OpenInventorXtExtendedViewer::AddButton (
 Widget aMenu
,const G4String& aLabel
,XtCallbackProc aCallback
)
{
  Widget widget = XmCreatePushButton(aMenu,(char*)aLabel.c_str(),NULL,0);
  XtManageChild(widget);
  XtAddCallback(widget,XmNactivateCallback,aCallback,(XtPointer)this);
}

void G4OpenInventorXtExtendedViewer::HelpCancelCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  XtUnmanageChild(This->fHelpForm);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void G4OpenInventorXtExtendedViewer::EscapeCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->Escape();
}

// Allow escape from X event loop via key
void G4OpenInventorXtExtendedViewer::EscapeFromKeyboardCbk(void* o) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)o;
  This->Escape();
}

void G4OpenInventorXtExtendedViewer::PostScriptCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  // FWJ Workaround: avoids empty 2nd page in file
  SbBool superimpState =
     This->fViewer->getSuperimpositionEnabled(This->fViewer->superimposition);
  This->fViewer->setSuperimpositionEnabled(This->fViewer->superimposition,
                                           FALSE);
  This->WritePostScript();
  if (superimpState)
     This->fViewer->setSuperimpositionEnabled(This->fViewer->superimposition,
                                              TRUE);
}
void G4OpenInventorXtExtendedViewer::PDFCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  // FWJ Workaround: avoids empty 2nd page in file
  SbBool superimpState =
     This->fViewer->getSuperimpositionEnabled(This->fViewer->superimposition);
  This->fViewer->setSuperimpositionEnabled(This->fViewer->superimposition,
                                           FALSE);
  This->WritePDF();
  if (superimpState)
     This->fViewer->setSuperimpositionEnabled(This->fViewer->superimposition,
                                              TRUE);
}

void G4OpenInventorXtExtendedViewer::PixmapPostScriptCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->WritePixmapPostScript();
}

void G4OpenInventorXtExtendedViewer::SceneGraphStatisticsCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SceneGraphStatistics();
}

void G4OpenInventorXtExtendedViewer::WriteInventorCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->WriteInventor();
}

void G4OpenInventorXtExtendedViewer::EraseDetectorCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->EraseDetector();
}

void G4OpenInventorXtExtendedViewer::EraseEventCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->EraseEvent();
}

void G4OpenInventorXtExtendedViewer::SetSolidCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SetSolid();
}

void G4OpenInventorXtExtendedViewer::SetWireFrameCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SetWireFrame();
}

void G4OpenInventorXtExtendedViewer::SetReducedWireFrameCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SetReducedWireFrame(true);
}

void G4OpenInventorXtExtendedViewer::SetFullWireFrameCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SetReducedWireFrame(false);
}

void G4OpenInventorXtExtendedViewer::UpdateSceneCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->UpdateScene();
}

void G4OpenInventorXtExtendedViewer::SetPreviewCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SetPreview();
}

void G4OpenInventorXtExtendedViewer::SetPreviewAndFullCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  This->SetPreviewAndFull();
}

void G4OpenInventorXtExtendedViewer::HelpCbk(
  Widget,XtPointer aData,XtPointer) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)aData;
  XtManageChild(This->fHelpForm);
  XmTextSetString(This->fHelpText,(char*)This->Help().c_str());
}
