//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VRML2Viewer.cc,v 1.9 2002-06-23 03:31:51 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2Viewer.cc
// Satoshi Tanaka & Yasuhide Sawada

#ifndef WIN32

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_VIEW

#include <math.h>

#include "G4Scene.hh"
#include "G4VRML2Viewer.hh"
#include "G4VRML2SceneHandler.hh"
#include "G4VRML2.hh"
#include "G4ios.hh"

G4VRML2Viewer::G4VRML2Viewer(G4VRML2SceneHandler& scene, const G4String& name) :
 G4VViewer(scene, scene.IncrementViewCount(), name),
 fSceneHandler(scene),
 fDest(scene.fDest)
{
	fViewHalfAngle = 0.5 * 0.785398 ; // 0.5 * 45*deg
	fsin_VHA       = sin ( fViewHalfAngle ) ;	
}

G4VRML2Viewer::~G4VRML2Viewer()
{}

void G4VRML2Viewer::SetView()
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML2Viewer::SetView(): No effects" << G4endl;
#endif

// Do nothing, since VRML a browser is running as a different process.
// SendViewParameters () will do this job instead.

}


void G4VRML2Viewer::DrawView()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4VRML2Viewer::DrawView()" << G4endl;
#endif
	fSceneHandler.VRMLBeginModeling() ;

        // Viewpoint node
        SendViewParameters(); 

	// Here is a minimal DrawView() function.
	NeedKernelVisit();
	ProcessView();
	FinishView();
}

void G4VRML2Viewer::ClearView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML2Viewer::ClearView(): No effects" << G4endl;
#endif
}

void G4VRML2Viewer::ShowView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML2Viewer::ShowView()" << G4endl;
#endif
	fSceneHandler.VRMLEndModeling();
}

void G4VRML2Viewer::FinishView(void)
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4VRML2Viewer::FinishView(): No effects" << G4endl;
#endif
}

void G4VRML2Viewer::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
      G4cerr << "***** G4VRML2Viewer::SendViewParameters()\n";
#endif 
	// error recovery
	if ( fsin_VHA < 1.0e-6 ) { return ; } 

	// camera distance
	G4double extent_radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
	G4double camera_distance = extent_radius / fsin_VHA ;

	// camera position on Z axis
	const G4Point3D&	target_point
	  = fSceneHandler.GetScene()->GetStandardTargetPoint()
	  + fVP.GetCurrentTargetPoint();
	G4double		E_z = target_point.z() + camera_distance;
	G4Point3D		E(0.0, 0.0, E_z );

	// VRML codes are generated below	
	fDest << "\n";
	fDest << "#---------- CAMERA" << "\n";
	fDest << "Viewpoint {"         << "\n";
	fDest << "\t" << "position "           ;
	fDest                 << E.x() << " "  ;
	fDest                 << E.y() << " "  ;
	fDest                 << E.z() << "\n" ;
	fDest << "}" << "\n";
	fDest << "\n";

} 



#endif
#endif //WIN32
