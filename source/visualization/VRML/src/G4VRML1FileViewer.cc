// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1FileViewer.cc,v 1.3 1999-11-04 02:38:48 stanaka Exp $
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
  G4cerr << "***** G4VRML1FileViewer::SetView(): No effects" << endl;
#endif
}

void G4VRML1FileViewer::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1FileViewer::DrawView()" << endl;
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
  G4cerr << "***** G4VRML1File1View::ClearView(): No effects" << endl;
#endif
}

void G4VRML1FileViewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1FileViewer::ShowView()" << endl;
#endif
	fSceneHandler.VRMLEndModeling();
}

void G4VRML1FileViewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML1FileViewer::FinishView(): No effects" << endl;
#endif
}


#endif
