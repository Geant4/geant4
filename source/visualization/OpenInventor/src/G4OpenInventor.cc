// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventor.cc,v 1.5 1999-05-12 14:00:58 barrand Exp $
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
#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorViewer.hh"


G4OpenInventor::G4OpenInventor (
 const G4String name
,const G4String nickname
,G4VGraphicsSystem::Functionality f
)
:G4VGraphicsSystem(name,nickname,f)
,interactorManager(NULL)
{
}

G4OpenInventor::~G4OpenInventor () {}

void G4OpenInventor::SetInteractorManager (G4VInteractorManager* im) {
  interactorManager = im;
}
G4VInteractorManager* G4OpenInventor::GetInteractorManager () {
  return interactorManager;
}
G4VSceneHandler* G4OpenInventor::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* p = new G4OpenInventorSceneHandler (*this, name);
  G4cout << G4OpenInventorSceneHandler::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    p;
}

G4VViewer* G4OpenInventor::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
  G4OpenInventorSceneHandler* pScene = (G4OpenInventorSceneHandler*)&scene;
  G4OpenInventorViewer*  pView  = new G4OpenInventorViewer (*pScene, name);
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
