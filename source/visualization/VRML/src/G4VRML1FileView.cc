// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1FileView.cc,v 1.1 1999-01-07 16:15:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRMLView.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include "G4SceneData.hh"
#include "G4VRML1FileView.hh"
#include "G4VRML1FileScene.hh"
#include "G4VRML1File.hh"
#include "G4ios.hh"

G4VRML1FileView::G4VRML1FileView(G4VRML1FileScene& scene,
				 const G4String& name) :
	G4VView(scene, scene.IncrementViewCount(), name), fScene(scene)
{}

G4VRML1FileView::~G4VRML1FileView()
{}

void G4VRML1FileView::SetView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1FileView::SetView()" << endl;
	G4cerr << "G4VRML1FileView::SetView(); not imlemented. " << endl;
#endif
}

void G4VRML1FileView::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1FileView::DrawView()" << endl;
#endif
	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML1FileView::ClearView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1File1View::ClearView()" << endl;
	G4cerr << "G4VRML1FileView::ClearView(); not implemented. " << endl;
#endif
}

void G4VRML1FileView::ShowView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1FileView::ShowView()" << endl;
#endif
	fScene.endSending();
}

void G4VRML1FileView::FinishView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML1FileView::FinishView()" << endl;
#endif
}


#endif
