// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1File.cc,v 1.1 1999-01-07 16:15:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1File.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


#include <stdio.h> // sscanf
#include <stdlib.h> // getenv

#include "G4VScene.hh"

#include "G4VRML1File.hh"
#include "G4VRML1FileScene.hh"
#include "G4VRML1FileView.hh"

#include "G4FRClient.hh"


G4VRML1File::G4VRML1File() :
	G4VGraphicsSystem("VRML1FILE", "VRML1FILE", G4VGraphicsSystem::threeD)
{
}

G4VRML1File::~G4VRML1File()
{
}


G4VScene* G4VRML1File::CreateScene(const G4String& name) 
{
	G4VScene *p = NULL;

	p = new G4VRML1FileScene(*this, name);

	G4cout << G4VRML1FileScene::GetSceneCount()
		<< " " << fName << " scenes extanct." << endl;

	return p;
}

G4VView* G4VRML1File::CreateView(G4VScene& scene, const G4String& name)
{
	G4VView* pView = NULL;

	G4VRML1FileScene* pScene = (G4VRML1FileScene*)&scene;
	pView = new G4VRML1FileView(*pScene, name);

	return pView;
}

#endif
