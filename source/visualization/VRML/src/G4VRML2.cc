// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2.cc,v 1.5 2000-08-19 18:34:59 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


#include <stdio.h> // sscanf
#include <stdlib.h> // getenv

#include "G4VSceneHandler.hh"

#include "G4VRML2.hh"
#include "G4VRML2SceneHandler.hh"
#include "G4VRML2Viewer.hh"

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


G4VSceneHandler* G4VRML2::CreateSceneHandler(const G4String& name) 
{
	G4VSceneHandler *p = NULL;

	p = new G4VRML2SceneHandler(*this, name);

	G4cout << G4VRML2SceneHandler::GetSceneCount()
		<< " " << fName << " scene handlers extanct." << G4endl;

	return p;
}

G4VViewer* G4VRML2::CreateViewer(G4VSceneHandler& scene, const G4String& name)
{
	G4VViewer* pView = NULL;

	G4VRML2SceneHandler* pScene = (G4VRML2SceneHandler*)&scene;
	pView = new G4VRML2Viewer(*pScene, name);

	return pView;
}

#endif
