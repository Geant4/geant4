// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1FileSceneHandler.cc,v 1.4 2000-04-27 13:56:16 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1FileSceneHandler.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


//#define DEBUG_FR_SCENE

#include "g4std/fstream"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "globals.hh"
#include "G4Scene.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Point3D.hh"
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

#include "G4VRML1FileSceneHandler.hh"
#include "G4VRML1FileViewer.hh"
#include "G4VRML1File.hh"

// CONST

const char  WRL_FILE_HEADER      [] = "g4_";
const char  DEFAULT_WRL_FILE_NAME[] = "g4.wrl";
const char  ENV_VRML_VIEWER      [] = "G4VRMLFILE_VIEWER";
const char  NO_VRML_VIEWER       [] = "NONE";
const char  VRMLFILE_DEST_DIR    [] = "G4VRMLFILE_DEST_DIR";


G4VRML1FileSceneHandler::G4VRML1FileSceneHandler(G4VRML1File& system, const G4String& name) :
	G4VSceneHandler(system, fSceneIdCount++, name),
	fSystem(system),
	fDest()   ,
	fFlagDestOpen( false ) 
{
	fSceneCount++;
	fCurrentDEF = "";
	strcpy(fVRMLFileName, "");

	if ( getenv( VRMLFILE_DEST_DIR ) == NULL ) {
		strcpy( fVRMLFileDestDir, "" );
	} else {
		strcpy( fVRMLFileDestDir, getenv( VRMLFILE_DEST_DIR ) );
	}

	// maximum number of g4.prim files in the dest directory
	fMaxFileNum = 1 ; // initialization
	if ( getenv( "G4VRMLFILE_MAX_FILE_NUM" ) != NULL ) {	
		
		sscanf( getenv("G4VRMLFILE_MAX_FILE_NUM"), "%d", &fMaxFileNum ) ;

	} else {
		fMaxFileNum = 1 ;
	}
	if( fMaxFileNum < 1 ) { fMaxFileNum = 1 ; }


}


G4VRML1FileSceneHandler::~G4VRML1FileSceneHandler()
{
#if defined DEBUG_FR_SCENE
	G4cerr << "***** ~G4VRML1FileSceneHandler" << G4endl;
#endif 
	fSceneCount--;
}


#define  G4VRML1SCENE   G4VRML1FileSceneHandler
#define  IS_CONNECTED   this->isConnected() 
#include "G4VRML1SceneHandlerFunc.icc"
#undef   IS_CONNECTED
#undef   G4VRML1SCENE 


void G4VRML1FileSceneHandler::connectPort()
{
	// g4.wrl, g4_1.wrl, ..., g4_MAX_FILE_INDEX.wrl
	const int MAX_FILE_INDEX = fMaxFileNum - 1 ;

	// dest directory (null if no environmental variables is set)
	strcpy ( fVRMLFileName, fVRMLFileDestDir) ; 

	// create full path name (default)
	strcat ( fVRMLFileName, DEFAULT_WRL_FILE_NAME );

	// Determine VRML file name
	for( int i = 0 ; i < fMaxFileNum ; i++) { 

		// Message
		if( fMaxFileNum > 1 && i == MAX_FILE_INDEX ) {
		  G4cerr << "==========================================="   << G4endl; 
		  G4cerr << "WARNING MESSAGE from VRML1FILE driver:      "  << G4endl;
		  G4cerr << "  This file name is the final one in the   "   << G4endl;
		  G4cerr << "  automatic updation of the output file name." << G4endl; 
		  G4cerr << "  You may overwrite an existing file of   "    << G4endl; 
                  G4cerr << "  the same name.                          "    << G4endl;
		  G4cerr << "==========================================="   << G4endl; 
		}

		// re-determine file to G4VRMLFILE_DEST_DIR/g4_i.wrl for i>0
		if( i >  0 ) { 
			sprintf( fVRMLFileName, "%s%s%d.wrl" , fVRMLFileDestDir,  WRL_FILE_HEADER, i );
		}

		// check validity of the file name
		G4std::ifstream  fin ; 
		fin.open(fVRMLFileName) ;
		if(!fin) { 
			// new file	
			fin.close();  // error recovery
			break; 
		} else { 
			// already exists (try next) 
			fin.close(); 
		} 

	} // for 

	// open a VRML 1.0 file with determined file name
	G4cerr << "==========================================="  << G4endl; 
	G4cerr << "Output VRML 1.0 file: " << fVRMLFileName << G4endl; 
	G4cerr << "Muximal number of file in the destination directory: " << fMaxFileNum << G4endl; 
	G4cerr << "  (Customizable as: setenv G4VRMLFILE_MAX_FILE_NUM number) " << G4endl;
	G4cerr << "===========================================" << G4endl; 

	fDest.open(fVRMLFileName) ;  fFlagDestOpen =  true ;
}


void G4VRML1FileSceneHandler::closePort()
{
	char command[256] ;
	char viewer [256] ; 
	strcpy( viewer, NO_VRML_VIEWER ); // initialization
	if( getenv( ENV_VRML_VIEWER ) ) {
		strcpy( viewer, getenv( ENV_VRML_VIEWER ) ) ;
	}

	// close VRML file	
	fDest.close();  fFlagDestOpen = false ;
	G4cerr << "*** VRML 1.0 File  " << fVRMLFileName << "  is generated." << G4endl;

	
	// Invoke viewer 

	if ( !strcmp(viewer, NO_VRML_VIEWER )) {
		G4cerr << "MESSAGE from VRML1FILE driver:"     << G4endl;
		G4cerr << "    Set an environmental variable  " ;
		G4cerr <<      ENV_VRML_VIEWER << G4endl;
		G4cerr << "    if you want to visualize the generated VRML file" << G4endl; 
		G4cerr << "    automatically.  For example, " << G4endl;
		G4cerr << "    setenv  " << ENV_VRML_VIEWER << "  vrweb " << G4endl;
	} else {
		sprintf( command, "%s %s", viewer, fVRMLFileName  );   
		system( command );
	}
}

G4int G4VRML1FileSceneHandler::fSceneIdCount = 0;
G4int G4VRML1FileSceneHandler::fSceneCount = 0;

#endif
