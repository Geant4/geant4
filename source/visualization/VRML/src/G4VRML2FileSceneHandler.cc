// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2FileSceneHandler.cc,v 1.1 1999-01-09 16:27:50 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2FileSceneHandler.cc
// Satoshi Tanaka & Yasuhide Sawada

//=================//
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
//=================//


//#define DEBUG_FR_SCENE

#include <fstream.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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

#include "G4VRML2FileSceneHandler.hh"
#include "G4VRML2FileViewer.hh"
#include "G4VRML2File.hh"

// CONST

const char  WRL_FILE_HEADER      [] = "g4_";
const char  DEFAULT_WRL_FILE_NAME[] = "g4.wrl";
const char  ENV_VRML_VIEWER      [] = "G4VRMLFILE_VIEWER";
const char  NO_VRML_VIEWER       [] = "NONE";
const char  VRMLFILE_DEST_DIR    [] = "G4VRMLFILE_DEST_DIR";


G4VRML2FileSceneHandler::G4VRML2FileSceneHandler(G4VRML2File& system, const G4String& name) :
	G4VSceneHandler(system, fSceneIdCount++, name),
	fSystem(system),
	fDest()   ,
	fFlagDestOpen( false ),
	fPVPickable  ( false )  
{
	fSceneCount++;

	// output file name
	strcpy(fVRMLFileName, "");

	// destination directory
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
	G4cerr << "***** ~G4VRML2FileSceneHandler" << endl;
#endif 
	fSceneCount--;
}


#define  G4VRML2SCENE   G4VRML2FileSceneHandler
#define  IS_CONNECTED   this->isConnected() 
#include "G4VRML2SceneHandlerFunc.icc"
#undef   IS_CONNECTED
#undef   G4VRML2SCENE 


void G4VRML2FileSceneHandler::connectPort()
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
		  G4cerr << "==========================================="   << endl; 
		  G4cerr << "WARNING MESSAGE from VRML2FILE driver:      "  << endl;
		  G4cerr << "  This file name is the final one in the   "   << endl;
		  G4cerr << "  automatic updation of the output file name." << endl; 
		  G4cerr << "  You may overwrite an existing file of   "    << endl; 
                  G4cerr << "  the same name.                          "    << endl;
		  G4cerr << "==========================================="   << endl; 
		}

		// re-determine file to G4VRMLFILE_DEST_DIR/g4_i.wrl for i>0
		if( i >  0 ) { 
			sprintf( fVRMLFileName, "%s%s%d.wrl" , fVRMLFileDestDir,  WRL_FILE_HEADER, i );
		}

		// check validity of the file name
		ifstream  fin ; 
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

	// open a VRML 2.0 file with determined file name
	G4cerr << "==========================================="  << endl; 
	G4cerr << "Output VRML 2.0 file: " <<    fVRMLFileName << endl; 
	G4cerr << "Muximal number of file in the destination directory: " << fMaxFileNum << endl; 
	G4cerr << "  (Customizable as: setenv G4VRMLFILE_MAX_FILE_NUM number) " << endl;
	G4cerr << "===========================================" << endl; 

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
	G4cerr << "*** VRML 2.0 File  " << fVRMLFileName << "  is generated." << endl;

	
	// Invoke viewer 

	if ( !strcmp(viewer, NO_VRML_VIEWER )) {
		G4cerr << "MESSAGE from VRML2FILE driver:"     << endl;
		G4cerr << "    Set an environmental variable  " ;
		G4cerr <<      ENV_VRML_VIEWER << endl;
		G4cerr << "    if you want to visualize the generated VRML file" << endl; 
		G4cerr << "    automatically.  For example, " << endl;
		G4cerr << "    setenv  " << ENV_VRML_VIEWER << "  vrwave " << endl;
	} else {
		sprintf( command, "%s %s", viewer, fVRMLFileName  );   
		system( command );
	}
}

G4int G4VRML2FileSceneHandler::fSceneIdCount = 0;
G4int G4VRML2FileSceneHandler::fSceneCount = 0;

#endif
