// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1SceneHandler.cc,v 1.2 1999-12-15 14:54:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1SceneHandler.cc
// Yasuhide Sawada and Satoshi Tanaka

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_SCENE

#include <unistd.h>
#include "g4std/fstream"

#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Point3D.hh"
#include "G4VisAttributes.hh"
#include "G4Transform.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyline.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"

#include "G4VRML1SceneHandler.hh"
#include "G4VRML1Viewer.hh"
#include "G4VRML1.hh"



G4VRML1SceneHandler::G4VRML1SceneHandler(G4VRML1& system, const G4String& name) :
	G4VSceneHandler(system, fSceneIdCount++, name),
	fSystem(system),
	fDest()
{
	fSceneCount++;
	fCurrentDEF = "";
}


G4VRML1SceneHandler::~G4VRML1SceneHandler()
{
#if defined DEBUG_FR_SCENE
	G4cerr << "***** ~G4VRML1SceneHandler" << G4endl;
#endif 
	fSceneCount--;
}



#define  G4VRML1SCENE  G4VRML1SceneHandler
#define  IS_CONNECTED  fDest.isConnected() 
#include "G4VRML1SceneHandlerFunc.icc"
#undef   IS_CONNECTED
#undef   G4VRML1SCENE 


void G4VRML1SceneHandler::connectPort(G4int max_trial)
{
	G4int trial = 0 ;
	int port = fSystem.getPort();
	for (trial = 0; !fDest.isConnected()&& trial < max_trial; trial++, port++ ) {
		if (fDest.connect( (const char * )fSystem.getHostName(), port)) {
		    // INET domain connection is established
			G4cerr << "*** GEANT4 is connected to port  ";
			G4cerr << fDest.getPort(); 
			G4cerr << " of server  " << fSystem.getHostName() << G4endl;
			break; 
		} else { 
			// Connection failed. Try the next port.
			G4cerr << "*** GEANT4 incremented targeting port to ";
			G4cerr << port << G4endl;
		}

		sleep (1);

	} // for

	if (!fDest.isConnected()) {
		G4cerr << "*** INET Connection failed. " << G4endl;
		G4cerr << "    Maybe, you have not invoked viewer  g4vrmlview  yet, " << G4endl;
		G4cerr << "    or too many viewers are already running in the " << G4endl;
		G4cerr << "    server host(" << fSystem.getHostName() << "). " << G4endl;
	}
}

void G4VRML1SceneHandler::closePort()
{
	fDest.close();
	G4cerr << "*** INET Connection closed. " << G4endl;
}


G4int G4VRML1SceneHandler::fSceneIdCount = 0;
G4int G4VRML1SceneHandler::fSceneCount = 0;

#endif
