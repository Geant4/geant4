// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Wo.cc,v 2.1 1998/11/06 13:42:08 allison Exp $
// GEANT4 tag $Name: geant4-00 $
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
#include "G4WoView.hh"
#include "G4GoScene.hh"
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
G4VScene* G4Wo::CreateScene (
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoScene* pScene = new G4GoScene (*this, name);
  return     pScene;
}
/***************************************************************************/
G4VView* G4Wo::CreateView (
 G4VScene& scene,
 const G4String& name
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4GoScene* pScene = (G4GoScene*)&scene;
  G4VView*   pView  = new G4WoView (*pScene, name);
  return     pView;
}


#endif

