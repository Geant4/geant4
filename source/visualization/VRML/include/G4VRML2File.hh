// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2File.hh,v 1.6 1999-12-15 14:54:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2File.hh
// Satoshi Tanaka & Yasuhide Sawada

#if defined (G4VIS_BUILD_VRMLFILE_DRIVER) || defined (G4VIS_USE_VRMLFILE)

#ifndef G4VRML2FILE_HH
#define G4VRML2FILE_HH

#include "G4VGraphicsSystem.hh"

class G4VSceneHandler;

class G4VRML2File: public G4VGraphicsSystem {
public:
	G4VRML2File(); 
	virtual ~G4VRML2File();
	G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
	G4VViewer*  CreateViewer(G4VSceneHandler&, const G4String& name = "");

};

#endif //G4VRML2File_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER

