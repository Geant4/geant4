// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Xo.cc,v 2.2 1998/11/06 13:42:10 allison Exp $
// GEANT4 tag $Name: geant4-00 $
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
#include "G4XoView.hh"
#include "G4GoScene.hh"
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
  XDisplayPutFileInResourceDatabase (XtDisplay(top),"G4Xo.xrm");
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
G4VScene* G4Xo::CreateScene (
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoScene* pScene = new G4GoScene (*this, name);
  return     pScene;
}
/***************************************************************************/
G4VView* G4Xo::CreateView (
 G4VScene& scene,
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoScene* pScene = (G4GoScene*)&scene;
  G4VView*   pView  = new G4XoView  (*pScene, name);
  return     pView;
}

#endif
