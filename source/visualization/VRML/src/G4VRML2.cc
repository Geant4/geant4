// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2.cc,v 1.1 1999-01-07 16:15:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


#include <stdio.h> // sscanf
#include <stdlib.h> // getenv

#include "G4VScene.hh"

#include "G4VRML2.hh"
#include "G4VRML2Scene.hh"
#include "G4VRML2View.hh"

#include "G4FRClient.hh"


G4VRML2::G4VRML2() :
	G4VGraphicsSystem("VRML2", "VRML2", G4VGraphicsSystem::threeD)
{
	// port number
	fPort = FR_VRML_DEFAULT_PORT;
	char *pport = getenv(FR_VRML_PORT_ENV);
	if (pport) {
		sscanf(pport, "%d", &fPort);
	}

	// host name
	fHostName = "localhost" ; // G4String::operator = ( const char* cs )
	char *phostname =  getenv(FR_VRML_HOST_NAME_ENV); 
	if (phostname) {
		fHostName = phostname;
	}
}

G4VRML2::~G4VRML2()
{
}


G4VScene* G4VRML2::CreateScene(const G4String& name) 
{
	G4VScene *p = NULL;

	p = new G4VRML2Scene(*this, name);

	G4cout << G4VRML2Scene::GetSceneCount()
		<< " " << fName << " scenes extanct." << endl;

	return p;
}

G4VView* G4VRML2::CreateView(G4VScene& scene, const G4String& name)
{
	G4VView* pView = NULL;

	G4VRML2Scene* pScene = (G4VRML2Scene*)&scene;
	pView = new G4VRML2View(*pScene, name);

	return pView;
}

#endif
