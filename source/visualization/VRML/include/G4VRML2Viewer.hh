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
// $Id: G4VRML2Viewer.hh,v 1.7 2002-06-23 03:31:43 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2Viewer.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifndef WIN32

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML2_VIEWER_HH
#define G4VRML2_VIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"

#include "G4FRClient.hh"

class G4VRML2SceneHandler;

class G4VRML2Viewer: public G4VViewer {
public:
	G4VRML2Viewer(G4VRML2SceneHandler& scene, const G4String& name = "");
	virtual ~G4VRML2Viewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView(); // Do nothing. SendViewParameters will do its job.
	void SendViewParameters ()  ;

private:
	G4VRML2SceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
	G4FRClient&   fDest ;	

	G4double      fViewHalfAngle ;	
	G4double      fsin_VHA       ;	
};

#endif //G4VRML2_VIEW_HH
#endif //G4VIS_BUILD_VRML_DRIVER
#endif //WIN32
