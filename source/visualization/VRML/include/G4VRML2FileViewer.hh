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
// $Id: G4VRML2FileViewer.hh,v 1.6 2001-07-27 22:33:13 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2FileViewer.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifndef G4VRML2FILE_VIEWER_HH
#define G4VRML2FILE_VIEWER_HH

#include "g4std/fstream"
#include "G4VViewer.hh"
#include "globals.hh"

class G4VRML2FileSceneHandler;

class G4VRML2FileViewer: public G4VViewer {
public:
	G4VRML2FileViewer(G4VRML2FileSceneHandler& scene, const G4String& name = "");
	virtual ~G4VRML2FileViewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView(); // Do nothing. SendViewParameters will do its job.
	void SendViewParameters ()  ;

private:
	G4VRML2FileSceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
	G4std::ofstream&         fDest ;

	G4double      fViewHalfAngle ;	
	G4double      fsin_VHA       ;	

};

#endif //G4VRML2File_VIEW_HH
