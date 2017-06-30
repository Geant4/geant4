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
// $Id: G4VRML2FileSceneHandler.cc 104289 2017-05-23 13:24:09Z gcosmo $
//
// G4VRML2FileSceneHandler.cc
// Satoshi Tanaka & Yasuhide Sawada


//#define DEBUG_FR_SCENE

#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VisManager.hh"
#include "G4Point3D.hh"
#include "G4VisAttributes.hh"
#include "G4VModel.hh"
#include "G4Scene.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyline.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"

#include "G4VRML2FileSceneHandler.hh"
#include "G4VRML2FileViewer.hh"
#include "G4VRML2File.hh"

// CONST

const char  WRL_FILE_HEADER      [] = "g4_";
const char  DEFAULT_WRL_FILE_NAME[] = "g4.wrl";
const char  ENV_VRML_VIEWER      [] = "G4VRMLFILE_VIEWER";
const char  NO_VRML_VIEWER       [] = "NONE";
const char  VRMLFILE_DEST_DIR    [] = "G4VRMLFILE_DEST_DIR";
const int   DEFAULT_MAX_WRL_FILE_NUM = 100 ;


G4VRML2FileSceneHandler::G4VRML2FileSceneHandler(G4VRML2File& system, const G4String& name) :
	G4VSceneHandler(system, fSceneIdCount++, name),
	fSystem(system),
	fFlagDestOpen( false ),
	fPVPickable  ( false ),
        fDest()
{
	// output file name
	strcpy(fVRMLFileName, "");

	// destination directory
	if ( getenv( VRMLFILE_DEST_DIR ) == NULL ) {
		strcpy( fVRMLFileDestDir, "" );
	} else {
		strcpy( fVRMLFileDestDir, getenv( VRMLFILE_DEST_DIR ) );
	}


	// maximum number of g4.prim files in the dest directory
	fMaxFileNum = DEFAULT_MAX_WRL_FILE_NUM ; // initialization
	if ( getenv( "G4VRMLFILE_MAX_FILE_NUM" ) != NULL ) {	
		
		sscanf( getenv("G4VRMLFILE_MAX_FILE_NUM"), "%d", &fMaxFileNum ) ;

	} else {
		fMaxFileNum = DEFAULT_MAX_WRL_FILE_NUM ;
	}
	if( fMaxFileNum < 1 ) { fMaxFileNum = 1; }


	// PV name pickability 	
	if( getenv( "G4VRML_PV_PICKABLE" ) != NULL ) {

		int is_pickable ;
		sscanf( getenv("G4VRML_PV_PICKABLE"), "%d", &is_pickable ) ;

		if ( is_pickable ) { SetPVPickability ( true ) ; }
	} 

	// PV Transparency
	SetPVTransparency ();

}


G4VRML2FileSceneHandler::~G4VRML2FileSceneHandler()
{
#if defined DEBUG_FR_SCENE
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** ~G4VRML2FileSceneHandler" << G4endl;
#endif 
	VRMLEndModeling();
}


#define  G4VRML2SCENEHANDLER   G4VRML2FileSceneHandler
#define  IS_CONNECTED   this->isConnected() 
#include "G4VRML2SceneHandlerFunc.icc"
#undef   IS_CONNECTED
#undef   G4VRML2SCENEHANDLER


void G4VRML2FileSceneHandler::connectPort()
{
	// g4_00.wrl, g4_01.wrl, ..., g4_MAX_FILE_INDEX.wrl
	const int MAX_FILE_INDEX = fMaxFileNum - 1 ;

	// dest directory (null if no environmental variables is set)
	strcpy ( fVRMLFileName, fVRMLFileDestDir) ; 

	// create full path name (default)
	strcat ( fVRMLFileName, DEFAULT_WRL_FILE_NAME );

	// Determine VRML file name
	for( int i = 0 ; i < fMaxFileNum ; i++) { 

		// Message in the final execution
		if( i == MAX_FILE_INDEX ) 
		{
		  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
		    G4cout << "==========================================="   << G4endl; 
		    G4cout << "WARNING MESSAGE from VRML2FILE driver:     "   << G4endl;
		    G4cout << "  This file name is the final one in the   "   << G4endl;
		    G4cout << "  automatic updation of the output file name." << G4endl; 
		    G4cout << "  You may overwrite existing files, i.e.   "   << G4endl; 
		    G4cout << "  g4_XX.wrl.                               "   << G4endl;
		    G4cout << "==========================================="   << G4endl; 
		  }
		}

		// re-determine file name as G4VRMLFILE_DEST_DIR/g4_XX.wrl 
		std::ostringstream filename;
		filename
		<< fVRMLFileDestDir << WRL_FILE_HEADER
		<< std::setw(2) << std::setfill('0') << i << ".wrl";
		strncpy(fVRMLFileName,filename.str().c_str(),sizeof(fVRMLFileName));

		// check validity of the file name
		std::ifstream  fin ; 
		fin.open(fVRMLFileName) ;
		if(!fin) { 
			// new file	
			fin.close();  
			break; 
		} else { 
			// already exists (try next) 
			fin.close(); 
		} 

	} // for 

	// open a VRML 2.0 file with determined file name
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	  G4cout << "==========================================="  << G4endl; 
	  G4cout << "Output VRML 2.0 file: " <<    fVRMLFileName << G4endl; 
	  G4cout << "Maximum number of files in the destination directory: " << fMaxFileNum << G4endl; 
	  G4cout << "  (Customizable with the environment variable: G4VRMLFILE_MAX_FILE_NUM) " << G4endl;
	  G4cout << "===========================================" << G4endl; 
	}
	fDest.open(fVRMLFileName) ;  fFlagDestOpen =  true ;
}


void G4VRML2FileSceneHandler::closePort()
{
	char command[256] ;
	char viewer [256] ; 
	strcpy( viewer, NO_VRML_VIEWER ); // initialization
	if( getenv( ENV_VRML_VIEWER ) ) {
		strcpy( viewer, getenv( ENV_VRML_VIEWER ) ) ;
	}

	// close VRML file	
	fDest.close();  fFlagDestOpen = false ;
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	      G4cout << "*** VRML 2.0 File  " << fVRMLFileName << "  is generated." << G4endl;

	
	// Invoke viewer 

	if ( !strcmp(viewer, NO_VRML_VIEWER )) {
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
		G4cout << "MESSAGE from VRML2FILE driver:"     << G4endl;
		G4cout << "    Set an environmental variable  " ;
		G4cout <<      ENV_VRML_VIEWER << G4endl;
		G4cout << "    if you want to visualize the generated VRML file" << G4endl; 
		G4cout << "    automatically.  For example, " << G4endl;
		G4cout << "    setenv  " << ENV_VRML_VIEWER << "  vrwave " << G4endl;
	  }
	} else {
		std::ostringstream ossCommand;
		ossCommand << viewer << ' ' << fVRMLFileName;
		strncpy(command,ossCommand.str().c_str(),sizeof(command));
		(void) system( command );
	}
}

G4int G4VRML2FileSceneHandler::fSceneIdCount = 0;
