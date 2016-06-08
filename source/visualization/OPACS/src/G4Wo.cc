// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Wo.cc,v 1.3 1999/01/11 00:47:32 allison Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// 
// Guy Barrand 04 November 1996
// Wo graphics system factory.

#ifdef G4VIS_BUILD_OPACS_DRIVER

//Co
#include <CPrinter.h>
#include <Wo.h>
//G4
#include "G4Xt.hh"
#include "G4WoViewer.hh"
#include "G4GoSceneHandler.hh"
//This
#include "G4Wo.hh"

static G4VInteractorManager* interactorManager = NULL;
/***************************************************************************/
G4Wo::G4Wo (
)
:G4VGraphicsSystem ("Wo",G4VGraphicsSystem::threeD)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*.........................................................................*/
{
  interactorManager = G4Xt::getInstance ();
}
/***************************************************************************/
G4Wo::~G4Wo (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
G4VInteractorManager* G4Wo::GetInteractorManager (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return interactorManager;
}
/***************************************************************************/
G4VSceneHandler* G4Wo::CreateSceneHandler (
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoSceneHandler* pScene = new G4GoSceneHandler (*this, name);
  return     pScene;
}
/***************************************************************************/
G4VViewer* G4Wo::CreateViewer (
 G4VSceneHandler& scene,
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoSceneHandler* pScene = (G4GoSceneHandler*)&scene;
  G4VViewer*   pView  = new G4WoViewer (*pScene, name);
  return     pView;
}


#endif

