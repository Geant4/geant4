// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2View.hh,v 1.1 1999-01-07 16:15:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2View.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML2_VIEW_HH
#define G4VRML2_VIEW_HH

#include "G4VView.hh"
#include "globals.hh"

#include "G4FRClient.hh"

class G4VRML2Scene;

class G4VRML2View: public G4VView {
public:
	G4VRML2View(G4VRML2Scene& scene, const G4String& name = "");
	~G4VRML2View();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView(); // Do nothing. SendViewParameters will do its job.
	void SendViewParameters ()  ;

private:
	G4VRML2Scene& fScene; // Reference to Graphics Scene for this view.
	G4FRClient&   fDest ;	

	G4double      fViewHalfAngle ;	
	G4double      fsin_VHA       ;	
};

#endif //G4VRML2_VIEW_HH
#endif //G4VIS_BUILD_VRML_DRIVER
