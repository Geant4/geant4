// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventor.cc,v 1.1 1999-01-07 16:15:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifdef G4VIS_BUILD_OI_DRIVER

//#include <HEPVis/SoDB.h>
#include <HEPVis/nodes/SoG4Box.h>
#include <HEPVis/nodes/SoG4Tubs.h>
#include <HEPVis/nodes/SoG4Cons.h>
#include <HEPVis/nodes/SoG4Trd.h>
#include <HEPVis/nodes/SoG4Trap.h>
#include <HEPVis/nodekits/SoDetectorTreeKit.h>

#include "G4OpenInventor.hh"
#include "G4OpenInventorScene.hh"
#include "G4OpenInventorView.hh"


G4OpenInventor::G4OpenInventor (
 const G4String name
,const G4String nickname
,G4VGraphicsSystem::Functionality f
)
:G4VGraphicsSystem(name,nickname,f)
,interactorManager(NULL)
{
}
void G4OpenInventor::SetInteractorManager (G4VInteractorManager* im) {
  interactorManager = im;
}
G4VInteractorManager* G4OpenInventor::GetInteractorManager () {
  return interactorManager;
}
G4VScene* G4OpenInventor::CreateScene (const G4String& name) {
  G4VScene* p = new G4OpenInventorScene (*this, name);
  G4cout << G4OpenInventorScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    p;
}

G4VView* G4OpenInventor::CreateView (G4VScene& scene, const G4String& name) 
{
  G4OpenInventorScene* pScene = (G4OpenInventorScene*)&scene;
  G4OpenInventorView*  pView  = new G4OpenInventorView (*pScene, name);
  if (pView -> GetOIVisualFound ()) {
  }
  else {
    delete pView;
    pView  = NULL;
  }
  return   pView;
}

void G4OpenInventor::InitHEPVis()
{
  // The below is too much for most loaders :
  //  HEPVis_SoDB::init();  

  SoG4Box::initClass();
  SoG4Tubs::initClass();
  SoG4Cons::initClass();
  SoG4Trd::initClass();
  SoG4Trap::initClass();
  SoDetectorTreeKit::initClass();
}



#endif
