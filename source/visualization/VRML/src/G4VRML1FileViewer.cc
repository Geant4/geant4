// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1FileViewer.cc,v 1.4 1999-12-15 14:54:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRMLView.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include "G4Scene.hh"
#include "G4VRML1FileViewer.hh"
#include "G4VRML1FileSceneHandler.hh"
#include "G4VRML1File.hh"
#include "G4ios.hh"

G4VRML1FileViewer::G4VRML1FileViewer(G4VRML1FileSceneHandler& scene,
				 const G4String& name) :
  G4VViewer(scene, scene.IncrementViewCount(), name), fSceneHandler(scene)
{}

G4VRML1FileViewer::~G4VRML1FileViewer()
{}

void G4VRML1FileViewer::SetView()
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1FileViewer::SetView(): No effects" << G4endl;
#endif
}

void G4VRML1FileViewer::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1FileViewer::DrawView()" << G4endl;
#endif

	fSceneHandler.VRMLBeginModeling();

	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML1FileViewer::ClearView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1File1View::ClearView(): No effects" << G4endl;
#endif
}

void G4VRML1FileViewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1FileViewer::ShowView()" << G4endl;
#endif
	fSceneHandler.VRMLEndModeling();
}

void G4VRML1FileViewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1FileViewer::FinishView(): No effects" << G4endl;
#endif
}


#endif
