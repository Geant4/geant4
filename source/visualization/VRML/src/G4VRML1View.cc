// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1View.cc,v 1.1 1999-01-07 16:15:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRMLView.cc
// Yasuhide Sawada & Satoshi Tanaka

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include "G4SceneData.hh"
#include "G4VRML1View.hh"
#include "G4VRML1Scene.hh"
#include "G4VRML1.hh"
#include "G4ios.hh"

G4VRML1View::G4VRML1View(G4VRML1Scene& scene, const G4String& name) :
	G4VView(scene, scene.IncrementViewCount(), name), fScene(scene)
{}

G4VRML1View::~G4VRML1View()
{}

void G4VRML1View::SetView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1View::SetView()" << endl;
	G4cerr << "G4VRML1View::SetView(); not imlemented. " << endl;
#endif
}

void G4VRML1View::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1View::DrawView()" << endl;
#endif
	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML1View::ClearView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1View::ClearView()" << endl;
	G4cerr << "G4VRML1View::ClearView(); not implemented. " << endl;
#endif
}

void G4VRML1View::ShowView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1View::ShowView()" << endl;
#endif
	fScene.endSending();
}

void G4VRML1View::FinishView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1View::FinishView()" << endl;
	//G4cerr << "G4VRML1View::FinishView(); not implemented. " << endl;
#endif
	//fScene.endSending();
}

#endif
