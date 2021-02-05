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
// Frederick Jones TRIUMF 30 NOV 2017

#ifdef G4VIS_BUILD_OIQT_DRIVER

// this :
#include "G4OpenInventorQt.hh"

#include "G4SoQt.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorQtViewer.hh"

#ifndef G4GMAKE
#include "moc_G4OpenInventorQt.cpp"
#endif

G4OpenInventorQt::G4OpenInventorQt()
  : G4OpenInventor("OpenInventorQt", "OIQt", G4VGraphicsSystem::threeD),
    fInited(false)
{
}

void G4OpenInventorQt::Initialize()
{
  if(fInited) return; //Done

  // FWJ DEBUG
  //  G4cout << "G4OpenInventorQt: SETINTERACTORMANAGER " << G4SoQt::getInstance()
  //         << G4endl;

  SetInteractorManager(G4SoQt::getInstance());

  // FWJ DEBUG
  //  G4cout << "G4OpenInventorQt: GETINTERACTORMANAGER " << GetInteractorManager()
  //         << G4endl;

  // FWJ from G4OIXt: these have no counterpart in Qt
  //  GetInteractorManager() ->
  //    RemoveDispatcher((G4DispatchFunction)XtDispatchEvent);
  //  GetInteractorManager() ->
  //    AddDispatcher((G4DispatchFunction)SoXt::dispatchEvent);
  //  Widget top = (Widget)GetInteractorManager()->GetMainInteractor();

  // FWJ apparently not used in OpenGLQt:
  // if(getenv("XENVIRONMENT")==NULL) {
  //   XrmDatabase database = XrmGetDatabase(QtDisplay(top));
  //   if(database!=NULL) {
  //     XrmPutLineResource(&database,"*topShadowColor:white");
  //     XrmPutLineResource(&database,"*bottomShadowColor:black");
  //     XrmPutLineResource(&database,"*foreground:black");
  //     XrmPutLineResource(&database,"*background:lightgrey");
  //     XrmPutLineResource(&database,"*borderColor:lightgrey");
  //     XrmPutLineResource(&database,"*fontList:-*-helvetica-bold-r-*-*-*-120-*-*-*-*-iso8859-1");
  //     XrmPutLineResource(&database,"*help_popup.title:Help");
  //     XrmPutLineResource(&database,"*helpCancel.labelString:Cancel");
  //     XrmPutLineResource(&database,"*helpText.editMode:multi_line_edit");
  //     XrmPutLineResource(&database,"*helpText.columns:60");
  //     XrmPutLineResource(&database,"*helpText.rows:20");
  //     XrmPutLineResource(&database,"*helpText.background:white");
  //     XrmPutLineResource(&database,"*helpText.fontList:*courier*-r-*--14-*");
  //     XrmPutLineResource(&database,"*helpText.maxLength:8000");
  //   }
  // }

  // FWJ for now, create an independent main window

  // NOW should be done in G4SoQt [public G4VInteractorManager]
  //QWidget* mainWin = SoQt::init("Geant4");

  // Note that GetMainInteractor() returns G4Interactor [typedef void*]
  //QWidget* mainWin = (QWidget*)(GetInteractorManager()->GetMainInteractor());

  //  FWJ DEBUG
  //  G4cout << "OpenInventorQt: OBTAINED SoQt main window " << mainWin << G4endl;
  //  G4cout << "SoQt::getTopLevelWidget()" << SoQt::getTopLevelWidget() << G4endl;

  // In parent G4OpenInventor
  InitNodes();

  fInited = true;
}

G4OpenInventorQt::~G4OpenInventorQt()
{
}

G4VViewer* G4OpenInventorQt::CreateViewer(G4VSceneHandler& scene,
                                          const G4String& name)
{
  auto pView = new G4OpenInventorQtViewer(static_cast<G4OpenInventorSceneHandler&>(scene), name);

  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cerr <<
      "G4OpenInventorQt::CreateViewer: ERROR flagged by negative"
      " view id in G4OpenInventorQtViewer creation."
      "\n Destroying view and returning null pointer."
      << G4endl;
      delete pView;
      pView = 0;
    }
  }
  if (!pView) {
    G4cerr <<
    "G4OpenInventorQt::CreateViewer: ERROR: null pointer on new G4OpenInventorQtViewer."
    << G4endl;
  }

  Initialize();
  
  return pView;
}

#endif
