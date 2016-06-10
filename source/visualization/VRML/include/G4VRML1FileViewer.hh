//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VRML1FileViewer.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// G4VRML1FileViewer.hh
// Satoshi Tanaka and Yasuhide Sawada

#ifndef G4VRML1FILE_VIEWER_HH
#define G4VRML1FILE_VIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"

class G4VRML1FileSceneHandler;

class G4VRML1FileViewer: public G4VViewer {
public:
	G4VRML1FileViewer(G4VRML1FileSceneHandler& scene, const G4String& name = "");
	virtual ~G4VRML1FileViewer();
	void ClearView();
	void DrawView();
	void ShowView();
	void FinishView();
private:
	void SetView();

private:
	G4VRML1FileSceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
};

#endif //G4VRML1File_VIEW_HH
