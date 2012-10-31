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
// $Id: G4OpenInventorXtViewer.cc,v 1.26 2010-11-10 17:53:22 allison Exp $
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
// Frederick Jones and TJR  October 2012
// Extended viewer based on G4OpenInventorXt.hh
// Uses G4OpenInventorXtExaminerViewer.

#ifdef G4VIS_BUILD_OIX_DRIVER

// this :
#include "G4OpenInventorXtExtendedViewer.hh"

#include <Inventor/nodes/SoSelection.h>
#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/viewers/SoXtFlyViewer.h>

#include <X11/StringDefs.h>
#include <X11/Shell.h>

#include <Xm/Xm.h>
#include <Xm/PushB.h>
#include <Xm/Form.h>
#include <Xm/CascadeB.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>

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
:G4OpenInventorXtViewer (sceneHandler, name)
,fViewer(0)
{
  if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "Window name: " << fName << G4endl;
}

G4OpenInventorXtExtendedViewer::~G4OpenInventorXtExtendedViewer()
{}

void G4OpenInventorXtExtendedViewer::Initialise() {
  
  G4String wName = fName;
  
  Widget parent = (Widget)fInteractorManager->GetParentInteractor ();
  int width = 600;
  int height = 600;

  if(!parent) {  
    // Check if user has specified an X-Windows-type geometry string...
    char str[32];

    G4String sgeometry = fVP.GetXGeometryString();
    if(sgeometry.empty()) {
      G4cout << "ERROR: Geometry string \""
             << sgeometry	   
             << "\" is empty.  Using \"600x600\"."
             << G4endl;
      width = 600;
      height = 600;  
      sprintf(str,"%dx%d",width,height);
      sgeometry = str;
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
    AddButton(menu,"PS (gl2ps)",PostScriptCbk);
    AddButton(menu,"PS (pixmap)",PixmapPostScriptCbk);
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

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Allow escape from X event loop via key
void G4OpenInventorXtExtendedViewer::EscapeFromKeyboardCbk(void* o) {
  G4OpenInventorXtExtendedViewer* This = (G4OpenInventorXtExtendedViewer*)o;
  This->Escape();
}

#endif
