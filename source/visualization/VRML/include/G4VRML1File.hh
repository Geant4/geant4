// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1File.hh,v 1.2 1999-01-09 16:27:30 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1File.hh
// Satoshi Tanaka & Yasuhide Sawada

#if defined (G4VIS_BUILD_VRMLFILE_DRIVER) || (G4VIS_USE_VRMLFILE)

#ifndef G4VRML1FILE_HH
#define G4VRML1FILE_HH

#include "G4VGraphicsSystem.hh"

class G4VSceneHandler;

class G4VRML1File: public G4VGraphicsSystem {
public:
	G4VRML1File(); 
	~G4VRML1File();
	G4VSceneHandler* CreateScene(const G4String& name = "");
	G4VViewer*  CreateView(G4VSceneHandler&, const G4String& name = "");

};

#endif //G4VRML1File_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER

