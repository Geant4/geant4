// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2File.hh,v 1.1 1999-01-07 16:15:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2File.hh
// Satoshi Tanaka & Yasuhide Sawada

#if defined (G4VIS_BUILD_VRMLFILE_DRIVER) || (G4VIS_USE_VRMLFILE)

#ifndef G4VRML2FILE_HH
#define G4VRML2FILE_HH

#include "G4VGraphicsSystem.hh"

class G4VScene;

class G4VRML2File: public G4VGraphicsSystem {
public:
	G4VRML2File(); 
	~G4VRML2File();
	G4VScene* CreateScene(const G4String& name = "");
	G4VView*  CreateView(G4VScene&, const G4String& name = "");

};

#endif //G4VRML2File_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER

