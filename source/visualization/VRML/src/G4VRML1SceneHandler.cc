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
// $Id: G4VRML1SceneHandler.cc 66870 2013-01-14 23:38:59Z adotti $
//
// G4VRML1SceneHandler.cc
// Yasuhide Sawada and Satoshi Tanaka

#ifndef WIN32

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


//#define DEBUG_FR_SCENE

#include <unistd.h>
#include <fstream>

#include "globals.hh"
#include "G4VisManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Point3D.hh"
#include "G4Scene.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyline.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
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
	fCurrentDEF = "";
}


G4VRML1SceneHandler::~G4VRML1SceneHandler()
{
#if defined DEBUG_FR_SCENE
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** ~G4VRML1SceneHandler" << G4endl;
#endif 
}



#define  G4VRML1SCENEHANDLER  G4VRML1SceneHandler
#define  IS_CONNECTED  fDest.isConnected() 
#include "G4VRML1SceneHandlerFunc.icc"
#undef   IS_CONNECTED
#undef   G4VRML1SCENEHANDLER


void G4VRML1SceneHandler::connectPort(G4int max_trial)
{
	G4int trial = 0 ;
	int port = fSystem.getPort();
	for (trial = 0; !fDest.isConnected()&& trial < max_trial; trial++, port++ ) {
		if (fDest.connect( (const char * )fSystem.getHostName(), port)) {
		    // INET domain connection is established
		  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
			G4cout << "*** GEANT4 is connected to port  ";
			G4cout << fDest.getPort(); 
			G4cout << " of server  " << fSystem.getHostName() << G4endl;
		  }
			break; 
		} else { 
			// Connection failed. Try the next port.
		  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
			G4cout << "*** GEANT4 incremented targeting port to ";
			G4cout << port << G4endl;
		  }
		}

		sleep (1);

	} // for

	if (!fDest.isConnected()) {
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
		G4cout << "*** INET Connection failed. " << G4endl;
		G4cout << "    Maybe, you have not invoked viewer  g4vrmlview  yet, " << G4endl;
		G4cout << "    or too many viewers are already running in the " << G4endl;
		G4cout << "    server host(" << fSystem.getHostName() << "). " << G4endl;
	  }
	}
}

void G4VRML1SceneHandler::closePort()
{
	fDest.close();
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	      G4cout << "*** INET Connection closed. " << G4endl;
}


G4int G4VRML1SceneHandler::fSceneIdCount = 0;

#endif
#endif
