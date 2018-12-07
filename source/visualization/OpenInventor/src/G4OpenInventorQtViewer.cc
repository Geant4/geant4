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

// Frederick Jones TRIUMF 07 November 2017

#ifdef G4VIS_BUILD_OIQT_DRIVER

// this :
#include "G4OpenInventorQtViewer.hh"

#include "G4OpenInventorQtExaminerViewer.hh"

#include <Inventor/nodes/SoSelection.h>

#include <Inventor/Qt/SoQt.h>
// FWJ these are needed (why?) to use flags in SoQtExaminerViewer constr.
#include <Inventor/Qt/viewers/SoQtViewer.h>
#include <Inventor/Qt/viewers/SoQtFullViewer.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>

//#include <Xm/Xm.h>
//#include <Xm/PushB.h>
//#include <Xm/Form.h>
//#include <Xm/CascadeB.h>
//#include <Xm/RowColumn.h>
//#include <Xm/Text.h>

#include "HEPVis/actions/SoGL2PSAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"
#include "G4VisManager.hh"

G4OpenInventorQtViewer::G4OpenInventorQtViewer(
   G4OpenInventorSceneHandler& sceneHandler, const G4String& name)
   : G4OpenInventorViewer(sceneHandler, name)
   , fViewer(0)
{
   // FWJ fName is in G4VViewer parent of G4OpenInventorViewer
   if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
       G4cout << "Window name: " << fName << G4endl;
}


void G4OpenInventorQtViewer::Initialise()
{

  QWidget* parent = SoQt::getTopLevelWidget();

  G4cout << "G4OIQtViewer: Creating G4OIQtExaminerViewer with parent " << 
     parent << G4endl;

  //  fViewer = new SoQtExaminerViewer(parent, "Geant4", TRUE);
  fViewer = new G4OpenInventorQtExaminerViewer(parent, "Geant4", TRUE);

  //  G4String wName = fName;
  //
  //  QWidget parent = (QWidget)fInteractorManager->GetParentInteractor();
  int width = 600;
  int height = 600;

  // FWJ not sure what this is for
  //     fInteractorManager->AddShell(fShell);

  // FWJ or this:
  //   } else {
  //    char* str = fInteractorManager->GetCreationString();
  //    if(str!=0) wName = str;
  //    fViewer = new SoQtExaminerViewer(parent,wName.c_str(),TRUE);
  //  }

  fViewer->setSize(SbVec2s(width, height));

  // Have a GL2PS render action :
  const SbViewportRegion& vpRegion = fViewer->getViewportRegion();
  fGL2PSAction = new SoGL2PSAction(vpRegion);
  fViewer->setGLRenderAction(fGL2PSAction);

  // Else :
  G4cout << "G4OpenInventorQtViewer: setting scene graph " << 
     fSoSelection << G4endl;
  G4cout << "G4OpenInventorQtViewer: getNumChildren " << 
     fSoSelection->getNumChildren() << G4endl;
  fViewer->setSceneGraph(fSoSelection);
  fViewer->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_ADD);
  fViewer->viewAll();
  fViewer->saveHomePosition();
  fViewer->setTitle(fName);
  fViewer->show();

  // This SHOULD invoke the event loop:
  //  if(fShell) {
  QWidget* mainWin = SoQt::getTopLevelWidget();
  G4cout << "G4OIQtViewer: calling SoQt::show on mainWin = " << mainWin 
         << G4endl;
  SoQt::show(mainWin);
  fInteractorManager->FlushAndWaitExecution();
  //  }
  fInteractorManager->SetCreatedInteractor(fViewer->getWidget());
}

G4OpenInventorQtViewer::~G4OpenInventorQtViewer()
{
  //  if(fShell) fInteractorManager->RemoveShell(fShell);
  if(fViewer) {
    fViewer->setSceneGraph(0);
    //FIXME : SGI : the below "delete" block things.
    //FIXME : CoinXt : the below "delete" crashe in ~SoXtRenderArea.
    //FIXME : delete fViewer;
  }
  //  if(fShell) XtDestroyWidget(fShell);
}

void G4OpenInventorQtViewer::FinishView()
{
  if(!fViewer) return;
  fViewer->viewAll();
  fViewer->saveHomePosition();
}

void G4OpenInventorQtViewer::SetView()
{
  G4OpenInventorViewer::SetView();
  if(!fViewer) return;
  // Background.
  G4Colour b = fVP.GetBackgroundColour();
  fViewer->setBackgroundColor
    (SbColor((float)b.GetRed(),(float)b.GetGreen(),(float)b.GetBlue()));
}


void G4OpenInventorQtViewer::ViewerRender()
{
  if(!fViewer) return;
  fViewer->render();
}

SoCamera* G4OpenInventorQtViewer::GetCamera () {
  if(!fViewer) return 0;
  return fViewer->getCamera();
}


// FWJ need new implementation in SoQt for the following...
//Widget G4OpenInventorQtViewer::AddMenu(
// Widget aMenuBar
//,const G4String& aName
//,const G4String& aLabel
//)
//{
  // // Pulldown menu :
  // Widget menu = XmCreatePulldownMenu(aMenuBar,(char*)aName.c_str(),NULL,0);
  // // Cascade button :
  // Arg args[2];
  // XmString cps = 
  //   XmStringLtoRCreate((char*)aLabel.c_str(),(char*)XmSTRING_DEFAULT_CHARSET);
  // XtSetArg (args[0],XmNlabelString,cps);
  // XtSetArg (args[1],XmNsubMenuId,menu);
  // Widget widget = XmCreateCascadeButton(aMenuBar,(char*)aName.c_str(),args,2);
  // XmStringFree (cps);
  // XtManageChild(widget);
  // return menu;
//}

// void G4OpenInventorQtViewer::AddButton (
//  Widget aMenu
// ,const G4String& aLabel
// ,XtCallbackProc aCallback
// )
// {
  // Widget widget = XmCreatePushButton(aMenu,(char*)aLabel.c_str(),NULL,0);
  // XtManageChild(widget);
  // XtAddCallback(widget,XmNactivateCallback,aCallback,(XtPointer)this);
//}

//void G4OpenInventorQtViewer::HelpCancelCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // XtUnmanageChild(This->fHelpForm);
//}


// void G4OpenInventorQtViewer::EscapeCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->Escape();
//}

// void G4OpenInventorQtViewer::PostScriptCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->WritePostScript();
//}

//void G4OpenInventorQtViewer::PixmapPostScriptCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->WritePixmapPostScript();
//}

//void G4OpenInventorQtViewer::SceneGraphStatisticsCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->SceneGraphStatistics();
//}

//void G4OpenInventorQtViewer::WriteInventorCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->WriteInventor();
//}

//void G4OpenInventorQtViewer::EraseDetectorCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->EraseDetector();
//}

//void G4OpenInventorQtViewer::EraseEventCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->EraseEvent();
//}

//void G4OpenInventorQtViewer::SetSolidCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->SetSolid();
//}

//void G4OpenInventorQtViewer::SetWireFrameCbk(
  // Widget,XtPointer aData,XtPointer) {
  // G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
  // This->SetWireFrame();
//}

// void G4OpenInventorQtViewer::SetReducedWireFrameCbk(
//   Widget,XtPointer aData,XtPointer) {
//   G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
//   This->SetReducedWireFrame(true);
// }

// void G4OpenInventorQtViewer::SetFullWireFrameCbk(
//   Widget,XtPointer aData,XtPointer) {
//   G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
//   This->SetReducedWireFrame(false);
// }

// void G4OpenInventorQtViewer::UpdateSceneCbk(
//   Widget,XtPointer aData,XtPointer) {
//   G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
//   This->UpdateScene();
// }

// void G4OpenInventorQtViewer::SetPreviewCbk(
//   Widget,XtPointer aData,XtPointer) {
//   G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
//   This->SetPreview();
// }

// void G4OpenInventorQtViewer::SetPreviewAndFullCbk(
//   Widget,XtPointer aData,XtPointer) {
//   G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
//   This->SetPreviewAndFull();
// }

// void G4OpenInventorQtViewer::HelpCbk(
//   Widget,XtPointer aData,XtPointer) {
//   G4OpenInventorQtViewer* This = (G4OpenInventorQtViewer*)aData;
//   XtManageChild(This->fHelpForm);
//   XmTextSetString(This->fHelpText,(char*)This->Help().c_str());
// }

#endif
