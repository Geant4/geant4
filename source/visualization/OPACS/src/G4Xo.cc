// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Xo.cc,v 1.5 1999/12/15 14:54:03 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 04 November 1996
// Wo graphics system factory.

#ifdef G4VIS_BUILD_OPACS_DRIVER

//#define DEBUG

#include <stddef.h>

//Co
#include <CPrinter.h>
//Xx
#include <XWidget.h>
//G4
#include "G4Xt.hh"
#include "G4XoViewer.hh"
#include "G4GoSceneHandler.hh"
//This
#include "G4Xo.hh"

static G4VInteractorManager* interactorManager = NULL;
/***************************************************************************/
G4Xo::G4Xo (
)
:G4VGraphicsSystem ("Xo",G4VGraphicsSystem::threeD)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*.........................................................................*/
{
  interactorManager = G4Xt::getInstance ();
  Widget top = (Widget)interactorManager->GetMainInteractor ();
  if(top==NULL) {
    CWarn       ("G4Xo : Unable to init Xt.\n");
    return;
  }
  XWidgetSetTop                     (top);
  //XDisplayPutFileInResourceDatabase (XtDisplay(top),"G4Xo.xrm");
}
/***************************************************************************/
G4Xo::~G4Xo (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
G4VInteractorManager* G4Xo::GetInteractorManager (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return interactorManager;
}
/***************************************************************************/
G4VSceneHandler* G4Xo::CreateSceneHandler (
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoSceneHandler* pScene = new G4GoSceneHandler (*this, name);
  return     pScene;
}
/***************************************************************************/
G4VViewer* G4Xo::CreateViewer (
 G4VSceneHandler& scene,
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoSceneHandler* pScene = (G4GoSceneHandler*)&scene;
  G4VViewer*   pView  = new G4XoViewer  (*pScene, name);
  return     pView;
}

#endif
