// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2File.cc,v 1.1 1999-01-07 16:15:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2File.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


#include <stdio.h> // sscanf
#include <stdlib.h> // getenv

#include "G4VScene.hh"

#include "G4VRML2File.hh"
#include "G4VRML2FileScene.hh"
#include "G4VRML2FileView.hh"

#include "G4FRClient.hh"


G4VRML2File::G4VRML2File() :
	G4VGraphicsSystem("VRML2FILE", "VRML2FILE", G4VGraphicsSystem::threeD)
{
}

G4VRML2File::~G4VRML2File()
{
}


G4VScene* G4VRML2File::CreateScene(const G4String& name) 
{
	G4VScene *p = NULL;

	p = new G4VRML2FileScene(*this, name);

	G4cout << G4VRML2FileScene::GetSceneCount()
		<< " " << fName << " scenes extanct." << endl;

	return p;
}

G4VView* G4VRML2File::CreateView(G4VScene& scene, const G4String& name)
{
	G4VView* pView = NULL;

	G4VRML2FileScene* pScene = (G4VRML2FileScene*)&scene;
	pView = new G4VRML2FileView(*pScene, name);

	return pView;
}

#endif
