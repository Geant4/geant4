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

//#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QFont>

#include "HEPVis/actions/SoGL2PSAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIQt.hh"

#include "G4SoQt.hh"

#ifndef G4GMAKE
#include "moc_G4OpenInventorQtViewer.cpp"
#endif

G4OpenInventorQtViewer::G4OpenInventorQtViewer(
   G4OpenInventorSceneHandler& sceneHandler, const G4String& name)
   : G4OpenInventorViewer(sceneHandler, name), fViewer(0)
{
   // FWJ fName is in G4VViewer parent of G4OpenInventorViewer
   if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
     G4cout << "Window name: " << fName << G4endl;
}


void G4OpenInventorQtViewer::Initialise()
{

   QWidget* parent = SoQt::getTopLevelWidget();

   // FWJ DEBUG
   //  G4cout << "G4OIQtViewer: Creating G4OIQtExaminerViewer with parent " <<
   //     parent << G4endl;

   fViewer = new G4OpenInventorQtExaminerViewer(parent, fName, TRUE);

   auto UI = G4UImanager::GetUIpointer();
   auto uiQt = dynamic_cast<G4UIQt*>(UI->GetG4UIWindow());

   // Moved this to G4OpenInventorQtExaminerViewer::afterRealizeHook()
   ///////////////////////////////////////////////////////////
   //
   // This explicitly sets the TabWidget as parent before addTab():
   //   if (uiQt) uiQt->AddTabWidget(parent, QString(fName));
   ///////////////////////////////////////////////////////////

   // Simpler: calls addTab(), but causes viewer parts to show (temporarily)
   // in the "Useful tips" page !!
   //   if (uiQt) uiQt->AddViewerTab(parent, fName);
   // Leaves an empty viewer window frame hanging around:
   //   if (uiQt) uiQt->AddTabWidget(fViewer->getWidget(), QString(fName));

   //  G4String wName = fName;
   //
   //  QWidget parent = (QWidget)fInteractorManager->GetParentInteractor();

   int width = fVP.GetWindowSizeHintX();
   int height = fVP.GetWindowSizeHintY();

   // FWJ not sure what this is for
   //     fInteractorManager->AddShell(fShell);

   // FWJ or this:
   //   } else {
   //    char* str = fInteractorManager->GetCreationString();
   //    if(str!=0) wName = str;
   //    fViewer = new SoQtExaminerViewer(parent,wName.c_str(),TRUE);
   //  }

   fViewer->setSize(SbVec2s(width, height));

   // Add common menu items...

   //   QMenuBar* menubar = fViewer->getMenubar();
   QMenu* filemenu = fViewer->getFileMenu();
   QMenu* etcmenu = fViewer->getEtcMenu();
   QFont* font = fViewer->getFont();

   // File menu

   FileWritePS = new QAction("Write PostScript (gl2ps)", this);
   FileWritePS->setFont(*font);
   connect(FileWritePS, SIGNAL(triggered()), this,
           SLOT(FileWritePSCB()));
   filemenu->addAction(FileWritePS);

   FileWritePDF = new QAction("Write PDF (gl2ps)", this);
   FileWritePDF->setFont(*font);
   connect(FileWritePDF, SIGNAL(triggered()), this,
           SLOT(FileWritePDFCB()));
   filemenu->addAction(FileWritePDF);

   FileWriteIV = new QAction("Write IV", this);
   FileWriteIV->setFont(*font);
   connect(FileWriteIV, SIGNAL(triggered()), this,
           SLOT(FileWriteIVCB()));
   filemenu->addAction(FileWriteIV);

   FileEscape = new QAction("Escape", this);
   FileEscape->setFont(*font);
   connect(FileEscape, SIGNAL(triggered()), this,
           SLOT(FileEscapeCB()));
   filemenu->addAction(FileEscape);

   //   G4cout << "G4OIQtViewer: externalApp = " << 
   //   static_cast<G4SoQt*>(fInteractorManager)->IsExternalApp() << G4endl;
   if (static_cast<G4SoQt*>(fInteractorManager)->IsExternalApp())
      fViewer->setExternalQtApp();

   // Register escape CB with viewer, allowing E key escape
   //   fViewer->addEscapeCallback(FileEscapeCB);
   //   fViewer->addEscapeCallback(FileEscapeCB, (void*)this);

   // Etc menu

   EtcEraseDetector = new QAction("Erase detector", this);
   EtcEraseDetector->setFont(*font);
   connect(EtcEraseDetector, SIGNAL(triggered()), this,
           SLOT(EtcEraseDetectorCB()));
   etcmenu->addAction(EtcEraseDetector);

   EtcEraseEvent = new QAction("Erase event", this);
   EtcEraseEvent->setFont(*font);
   connect(EtcEraseEvent, SIGNAL(triggered()), this,
           SLOT(EtcEraseEventCB()));
   etcmenu->addAction(EtcEraseEvent);

   EtcSetSolid = new QAction("Set solid", this);
   EtcSetSolid->setFont(*font);
   connect(EtcSetSolid, SIGNAL(triggered()), this, SLOT(EtcSetSolidCB()));
   etcmenu->addAction(EtcSetSolid);

   EtcSetReducedWireframe = new QAction("Set (G4) reduced wireframe", this);
   EtcSetReducedWireframe->setFont(*font);
   connect(EtcSetReducedWireframe, SIGNAL(triggered()), this,
           SLOT(EtcSetReducedWireframeCB()));
   etcmenu->addAction(EtcSetReducedWireframe);

   EtcSetFullWireframe = new QAction("Set full wireframe", this);
   EtcSetFullWireframe->setFont(*font);
   connect(EtcSetFullWireframe, SIGNAL(triggered()), this,
           SLOT(EtcSetFullWireframeCB()));
   etcmenu->addAction(EtcSetFullWireframe);

   EtcVisibMInvisibD = new QAction("Visible mothers + invisible daughters",
                                   this);
   EtcVisibMInvisibD->setFont(*font);
   connect(EtcVisibMInvisibD, SIGNAL(triggered()), this,
           SLOT(EtcVisibMInvisibDCB()));
   etcmenu->addAction(EtcVisibMInvisibD);

   EtcVisibMVisibD = new QAction("Visible mothers + visible daughters", this);
   EtcVisibMVisibD->setFont(*font);
   connect(EtcVisibMVisibD, SIGNAL(triggered()), this,
           SLOT(EtcVisibMVisibDCB()));
   etcmenu->addAction(EtcVisibMVisibD);

   EtcUpdateScene = new QAction("Update scene", this);
   EtcUpdateScene->setFont(*font);
   connect(EtcUpdateScene, SIGNAL(triggered()), this,
           SLOT(EtcUpdateSceneCB()));
   etcmenu->addAction(EtcUpdateScene);

   EtcSceneGraphStats = new QAction("Scene graph stats", this);
   EtcSceneGraphStats->setFont(*font);
   connect(EtcSceneGraphStats, SIGNAL(triggered()), this,
           SLOT(EtcSceneGraphStatsCB()));
   etcmenu->addAction(EtcSceneGraphStats);


  // Have a GL2PS render action :
  const SbViewportRegion& vpRegion = fViewer->getViewportRegion();
  fGL2PSAction = new SoGL2PSAction(vpRegion);
  fViewer->setGLRenderAction(fGL2PSAction);

  // Else :

  // FWJ DEBUG
  //  G4cout << "G4OpenInventorQtViewer: setting scene graph " <<
  //     fSoSelection << G4endl;
  //  G4cout << "G4OpenInventorQtViewer: getNumChildren " <<
  //     fSoSelection->getNumChildren() << G4endl;

  fViewer->setSceneGraph(fSoSelection);
  fViewer->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_ADD);
  fViewer->viewAll();
  fViewer->saveHomePosition();
  // SOMEHOW this also the OIQt main window title
  if (!uiQt) fViewer->setTitle(fName);
  fViewer->show();

  // This SHOULD invoke the event loop:
  //  if(fShell) {

  QWidget* mainWin = SoQt::getTopLevelWidget();

  // FWJ DEBUG
  //  G4cout << "G4OIQtViewer: calling SoQt::show on mainWin = " << mainWin
  //         << G4endl;

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


// File menu...

void G4OpenInventorQtViewer::FileWritePSCB()
{
   //   G4cout << "G4OIQtViewer: File: Write PS CALLBACK" << G4endl;
   // FWJ Workaround: avoids empty 2nd page in file
   SbBool superimpState =
      fViewer->getSuperimpositionEnabled(fViewer->superimposition);
   fViewer->setSuperimpositionEnabled(fViewer->superimposition, FALSE);
   WritePostScript();
   if (superimpState)
      fViewer->setSuperimpositionEnabled(fViewer->superimposition, TRUE);
}

void G4OpenInventorQtViewer::FileWritePDFCB()
{
   //   G4cout << "G4OIQtViewer: File: Write PDF CALLBACK" << G4endl;
   // FWJ Workaround: avoids empty 2nd page in file
   SbBool superimpState =
      fViewer->getSuperimpositionEnabled(fViewer->superimposition);
   fViewer->setSuperimpositionEnabled(fViewer->superimposition, FALSE);
   WritePDF();
   if (superimpState)
      fViewer->setSuperimpositionEnabled(fViewer->superimposition, TRUE);
}

void G4OpenInventorQtViewer::FileWriteIVCB()
{
   //   G4cout << "G4OIQtViewer: File: Write IV CALLBACK" << G4endl;
   WriteInventor();
}

void G4OpenInventorQtViewer::FileEscapeCB()
{
   //   G4cout << "G4OIQtViewer: File: Escape CALLBACK" << G4endl;
   static_cast<G4SoQt*>(fInteractorManager)->ExitSecondaryLoop();
   //   Escape();
}

// Etc menu...

void
G4OpenInventorQtViewer::EtcEraseDetectorCB()
{
   //   G4cout << "G4OIQtViewer: Etc: Erase Detector CALLBACK" << G4endl;
   EraseDetector();
}

void
G4OpenInventorQtViewer::EtcEraseEventCB()
{
   //   G4cout << "G4OIQtViewer: Etc: Erase Event CALLBACK" << G4endl;
   EraseEvent();
}

void G4OpenInventorQtViewer::EtcSetSolidCB()
{
   //   G4cout << "G4OIQtViewer: Etc: Set Solid CALLBACK" << G4endl;
   SetSolid();
}

void G4OpenInventorQtViewer::EtcSetReducedWireframeCB()
{
   // G4cout << "G4OIQtViewer: Etc: Set Reduced Wireframe CALLBACK" << G4endl;
   SetReducedWireFrame(true);
}

void G4OpenInventorQtViewer::EtcSetFullWireframeCB()
{
   // G4cout << "G4OIQtViewer: Etc: Set Full Wireframe CALLBACK" << G4endl;
   SetReducedWireFrame(false);
}

void G4OpenInventorQtViewer::EtcVisibMInvisibDCB()
{
   // G4cout << "G4OIQtViewer: Etc: Visible Mothers + Invisible Daughters"
   //   " CALLBACK" << G4endl;
   SetPreview();
}

void G4OpenInventorQtViewer::EtcVisibMVisibDCB()
{
   // G4cout << "G4OIQtViewer: Etc: Visible Mothers + Visible Daughters"
   // "CALLBACK" << G4endl;
   SetPreviewAndFull();
}

void G4OpenInventorQtViewer::EtcUpdateSceneCB()
{
   //   G4cout << "G4OIQtViewer: Etc: Update Scene CALLBACK" << G4endl;
   UpdateScene();
}

void G4OpenInventorQtViewer::EtcSceneGraphStatsCB()
{
   //   G4cout << "G4OIQtViewer: Etc: Scene Graph Stats CALLBACK" << G4endl;
   SceneGraphStatistics();
}


#endif
