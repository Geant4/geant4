// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1Viewer.cc,v 1.3 1999-11-04 02:38:48 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRMLView.cc
// Yasuhide Sawada & Satoshi Tanaka

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include "G4Scene.hh"
#include "G4VRML1Viewer.hh"
#include "G4VRML1SceneHandler.hh"
#include "G4VRML1.hh"
#include "G4ios.hh"

G4VRML1Viewer::G4VRML1Viewer(G4VRML1SceneHandler& scene, const G4String& name) :
   G4VViewer(scene, scene.IncrementViewCount(), name), fSceneHandler(scene)
{}

G4VRML1Viewer::~G4VRML1Viewer()
{}

void G4VRML1Viewer::SetView()
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1Viewer::SetView(): No effects" << endl;
#endif
}

void G4VRML1Viewer::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1Viewer::DrawView()" << endl;
#endif

	fSceneHandler.VRMLBeginModeling();

	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML1Viewer::ClearView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1Viewer::ClearView(): No effects" << endl;
#endif
}

void G4VRML1Viewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1Viewer::ShowView()" << endl;
#endif
	fSceneHandler.VRMLEndModeling();
}

void G4VRML1Viewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1Viewer::FinishView(): No effects" << endl;
#endif
}

#endif
