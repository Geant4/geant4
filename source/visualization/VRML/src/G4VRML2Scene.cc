// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2Scene.cc,v 1.1 1999-01-07 16:15:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2Scene.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_SCENE

#include <unistd.h>
#include <fstream.h>
#include <math.h>

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

#include "G4VRML2Scene.hh"
#include "G4VRML2View.hh"
#include "G4VRML2.hh"



G4VRML2Scene::G4VRML2Scene(G4VRML2& system, const G4String& name) :
	G4VScene(system, fSceneIdCount++, name),
	fSystem(system),
	fDest(), 
	fPVPickable ( false ) 
{
	fSceneCount++;

	
	// PV name pickability 	
	if( getenv( "G4VRML_PV_PICKABLE" ) != NULL ) {

		int is_pickable ;
		sscanf( getenv("G4VRML_PV_PICKABLE"), "%d", &is_pickable ) ;

		if ( is_pickable ) { SetPVPickability ( true ) ; }
	} 
	
	// PV Transparency
	SetPVTransparency ();

}


G4VRML2Scene::~G4VRML2Scene()
{
#if defined DEBUG_FR_SCENE
	G4cerr << "***** ~G4VRML2Scene" << endl;
#endif 
	fSceneCount--;
}



#define  G4VRML2SCENE  G4VRML2Scene
#define  IS_CONNECTED  fDest.isConnected() 
#include "G4VRML2SceneFunc.icc"
#undef   IS_CONNECTED
#undef   G4VRML2SCENE 


void G4VRML2Scene::connectPort(G4int max_trial)
{
	G4int trial = 0 ;
	int port = fSystem.getPort();
	for (trial = 0; !fDest.isConnected()&& trial < max_trial; trial++, port++ ) {
		if (fDest.connect( (const char * )fSystem.getHostName(), port)) {
		    // INET domain connection is established
			G4cerr << "*** GEANT4 is connected to port  ";
			G4cerr << fDest.getPort(); 
			G4cerr << " of server  " << fSystem.getHostName() << endl;
			break; 
		} else { 
			// Connection failed. Try the next port.
			G4cerr << "*** GEANT4 incremented targeting port to ";
			G4cerr << port << endl;
		}

		sleep (1);

	} // for

	if (!fDest.isConnected()) {
		G4cerr << "*** INET Connection failed. " << endl;
		G4cerr << "    Maybe, you have not invoked viewer  g4vrmlview  yet, " << endl;
		G4cerr << "    or too many viewers are already running in the " << endl;
		G4cerr << "    server host(" << fSystem.getHostName() << "). " << endl;
	}
}

void G4VRML2Scene::closePort()
{
	fDest.close();
	G4cerr << "*** INET Connection closed. " << endl;
}


G4int G4VRML2Scene::fSceneIdCount = 0;
G4int G4VRML2Scene::fSceneCount = 0;

#endif
