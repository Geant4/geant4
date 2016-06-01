// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayX.cc,v 2.2 1998/11/06 13:41:57 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing graphics system factory.

#ifdef G4VIS_BUILD_RAYX_DRIVER

#include "G4RayX.hh"

#include "G4VisFeaturesOfRay.hh"
#include "G4RayXMessenger.hh"
#include "G4RayScene.hh"
#include "G4RayXView.hh"

G4RayX::G4RayX ():
  G4VGraphicsSystem ("G4RayX",
		     "RayX",
		     G4VisFeaturesOfRayX (),
		     G4VGraphicsSystem::twoD) {
    if (fInstances >= 1) {
      G4Exception ("An attempt to instantiate a second G4RayX object.");
    }
    fInstances++;
    // Instantiate my messenger.
    fpMessenger = new G4RayXMessenger (this);
}

G4RayX::~G4RayX () {
  delete fpMessenger;
  fInstances--;
}

G4VScene* G4RayX::CreateScene (const G4String& name) {
  G4RayScene* pRayScene = new G4RayScene (*this, name);
  G4VScene* pScene = pRayScene;
  G4cout << G4RayScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  fpMessenger -> RegisterScene ((G4RayScene*) pScene);
  return pScene;
}

G4VView* G4RayX::CreateView (G4VScene& scene, const G4String& name) {
  G4RayView* pRayXView = new G4RayXView ((G4RayScene&) scene, name);
  G4VView* pView = pRayXView;
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "Error flagged by negative view id in G4RayXView creation."
	"\nDestroying view and returning null pointer."
	   <<endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "Null pointer on new G4RayXView." << endl;
  }
  fpMessenger -> RegisterView (pRayXView);
  return pView;
}

G4int G4RayX::fInstances = 0;

#endif
