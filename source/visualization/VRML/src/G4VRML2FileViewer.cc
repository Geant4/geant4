// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2FileViewer.cc,v 1.2 1999-01-11 00:48:11 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2FileViewer.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include <math.h>

#include "G4Scene.hh"
#include "G4VRML2FileViewer.hh"
#include "G4VRML2FileSceneHandler.hh"
#include "G4VRML2File.hh"
#include "G4ios.hh"

G4VRML2FileViewer::G4VRML2FileViewer(G4VRML2FileSceneHandler& scene,
				 const G4String& name) :
 G4VViewer(scene, scene.IncrementViewCount(), name),
 fSceneHandler(scene),
 fDest(scene.fDest)
{
	fViewHalfAngle = 0.5 * 0.785398 ; // 0.5 * 45*deg
	fsin_VHA       = sin ( fViewHalfAngle ) ;	
}

G4VRML2FileViewer::~G4VRML2FileViewer()
{}

void G4VRML2FileViewer::SetView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML2FileViewer::SetView()" << endl;
	G4cerr << "G4VRML2FileViewer::SetView(); not imlemented. " << endl;
#endif

// Do nothing, since VRML a browser is running as a different process.
// SendViewParameters () will do this job instead.

}

void G4VRML2FileViewer::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML2FileViewer::DrawView()" << endl;
#endif

	// Open VRML2 file and output header comments
	fSceneHandler.beginSending() ; 

        // Viewpoint node
        SendViewParameters(); 

	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML2FileViewer::ClearView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML2File1View::ClearView()" << endl;
	G4cerr << "G4VRML2FileViewer::ClearView(); not implemented. " << endl;
#endif
}

void G4VRML2FileViewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML2FileViewer::ShowView()" << endl;
#endif
	fSceneHandler.endSending();
}

void G4VRML2FileViewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML2FileViewer::FinishView()" << endl;
#endif
}

void G4VRML2FileViewer::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
      G4cerr << "***** G4VRML2FileViewer::SendViewParameters()\n";
#endif 

	// error recovery
	if ( fsin_VHA < 1.0e-6 ) { return ; } 

	// camera distance
	G4double extent_radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
	G4double camera_distance = extent_radius / fsin_VHA ;

	// camera position on Z axis
	const G4Point3D&	target_point = fVP.GetCurrentTargetPoint();
	G4double		E_z = target_point.z() + camera_distance;
	G4Point3D		E(0.0, 0.0, E_z );

	// VRML codes are generated below	
	fDest << endl;
	fDest << "#---------- CAMERA" << endl;
	fDest << "Viewpoint {"         << endl;
	fDest << "\t" << "position "           ;
	fDest                 << E.x() << " "  ;
	fDest                 << E.y() << " "  ;
	fDest                 << E.z() << endl ;
	fDest << "}" << endl;
	fDest << endl;

} 

#endif
