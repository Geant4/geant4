// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DAWNFILESceneHandler.cc,v 1.3 1999-12-15 14:54:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA
// DAWNFILE scene.


//=================//
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
//=================//


#define __G_ANSI_C__

// #define DEBUG_FR_SCENE

     //----- header files
#include "g4std/fstream"
#include <stdlib.h>
#include <string.h>
#include "globals.hh"
#include "G4FRConst.hh"
#include "G4DAWNFILE.hh"
#include "G4DAWNFILESceneHandler.hh"
#include "G4DAWNFILEViewer.hh"
#include "G4Point3D.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyline.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4ModelingParameters.hh"
#include "G4VPhysicalVolume.hh"

//----- constants
const char  FR_ENV_CULL_INVISIBLE_OBJECTS [] = "G4DAWN_CULL_INVISIBLE_OBJECTS";
const char  G4PRIM_FILE_HEADER      [] = "g4_";
const char  DEFAULT_G4PRIM_FILE_NAME[] = "g4.prim";


///////////////////////////
// Driver-dependent part //
///////////////////////////


	//----- G4DAWNFILESceneHandler, constructor
G4DAWNFILESceneHandler::G4DAWNFILESceneHandler (G4DAWNFILE& system, const G4String& name):
fSystem   (system)                        ,
G4VSceneHandler  (system, fSceneIdCount++, name) ,
fPrimDest ()                              ,
FRflag_in_modeling     (false)            ,
flag_saving_g4_prim    (false)            ,
COMMAND_BUF_SIZE       (G4FRofstream::SEND_BUFMAX)
{
	// count instantiated scenes
	fSceneCount++;

	// g4.prim filename and its directory
	if ( getenv( "G4DAWNFILE_DEST_DIR" ) == NULL ) {
		strcpy( fG4PrimDestDir , "" );
		strcpy( fG4PrimFileName, DEFAULT_G4PRIM_FILE_NAME );
	} else {
		strcpy( fG4PrimDestDir , getenv( "G4DAWNFILE_DEST_DIR" ) );
		strcpy( fG4PrimFileName, DEFAULT_G4PRIM_FILE_NAME        );
	}
		
	// maximum number of g4.prim files in the dest directory
	fMaxFileNum = 1 ; // initialization
	if ( getenv( "G4DAWNFILE_MAX_FILE_NUM" ) != NULL ) {	
		
		sscanf( getenv("G4DAWNFILE_MAX_FILE_NUM"), "%d", &fMaxFileNum ) ;

	} else {
		fMaxFileNum = 1 ;
	}
	if( fMaxFileNum < 1 ) { fMaxFileNum = 1 ; }

} 


	//----- G4DAWNFILESceneHandler, destructor
G4DAWNFILESceneHandler::~G4DAWNFILESceneHandler () 
{
#if defined DEBUG_FR_SCENE
	G4cerr << "***** ~G4DAWNFILESceneHandler" << G4endl;
#endif 
  fSceneCount--;
  ClearStore (); // clear current scene

}

//-----
void	G4DAWNFILESceneHandler::SetG4PrimFileName() 
{
	// g4.prim, g4_1.prim, ..., g4_MAX_FILE_INDEX.prim
	const int MAX_FILE_INDEX = fMaxFileNum - 1 ;

	// dest directory (null if no environmental variables is set)
	strcpy ( fG4PrimFileName, fG4PrimDestDir) ; 

	// create full path name (default)
	strcat ( fG4PrimFileName, DEFAULT_G4PRIM_FILE_NAME );

	// Automatic updation of file names
	for( int i = 0 ; i < fMaxFileNum ; i++) { 

		// Message
		if( fMaxFileNum > 1 && i == MAX_FILE_INDEX ) {
		  G4cerr << "==========================================="   << G4endl; 
		  G4cerr << "WARNING MESSAGE from DAWNFILE driver:      "   << G4endl;
		  G4cerr << "  This file name is the final one in the   "   << G4endl;
		  G4cerr << "  automatic updation of the output file name." << G4endl; 
		  G4cerr << "  You may overwrite an existing file of   "    << G4endl; 
                  G4cerr << "  the same name.                          "    << G4endl;
		  G4cerr << "==========================================="   << G4endl; 
		}

		// re-determine file to G4DAWNFILE_DEST_DIR/g4_i.prim for i>0
		if( i >  0 ) { 
			sprintf( fG4PrimFileName, "%s%s%d.prim" , fG4PrimDestDir,  G4PRIM_FILE_HEADER, i );
		}

		// check validity of the file name
		G4std::ifstream  fin ; 
		fin.open(fG4PrimFileName) ;
		if(!fin) { 
			// new file	
			fin.close();  // error recovery
			break; 
		} else { 
			// already exists (try next) 
			fin.close(); 
		} 

	} // for 

	G4cerr << "===========================================" << G4endl; 
	G4cerr << "Output file: " <<    fG4PrimFileName         << G4endl; 
	G4cerr << "Muximal number of file in the destination directory: " << fMaxFileNum << G4endl; 
	G4cerr << "  (Customizable as: setenv G4DAWNFILE_MAX_FILE_NUM number) " << G4endl;
	G4cerr << "===========================================" << G4endl; 

} // G4DAWNFILESceneHandler::SetG4PrimFileName()


//-----
void	G4DAWNFILESceneHandler::BeginSavingG4Prim( void ) 
{
#if defined DEBUG_FR_SCENE
	G4cerr << "***** BeginSavingG4Prim (called)\n";
#endif

	if( !IsSavingG4Prim() ) 
	{ 
#if defined DEBUG_FR_SCENE
	        G4cerr << "*****                   (started) " ;
	        G4cerr << "(open g4.prim, ##)"  << G4endl;
#endif
		SetG4PrimFileName() ; // returned to fG4PrimFileName
		fPrimDest.Open(fG4PrimFileName)   ;

		SendStr( FR_G4_PRIM_HEADER   )    ; 
		flag_saving_g4_prim = true        ; 
	} 
}

void	G4DAWNFILESceneHandler::EndSavingG4Prim  ( void ) 
{
#if defined DEBUG_FR_SCENE
	G4cerr << "***** EndSavingG4Prim (called)\n";
#endif

	if(  IsSavingG4Prim() )
	{
#if defined DEBUG_FR_SCENE
          G4cerr << "*****                 (started) (close g4.prim)" << G4endl;
#endif
		fPrimDest.Close()               ;
		flag_saving_g4_prim = false ; 
	} 
}


//----- 
void G4DAWNFILESceneHandler::FRBeginModeling( void )
{
	if( !FRIsInModeling() )  	
	{
#if defined DEBUG_FR_SCENE
	  G4cerr << "***** G4DAWNFILESceneHandler::FRBeginModeling (called & started)" << G4endl;
#endif

			//----- Send saving command and heading comment
		BeginSavingG4Prim();

			//----- Send bounding box command
		SendBoundingBox();

			//----- send SET_CAMERA command 
#if defined DEBUG_FR_SCENE
		G4cerr << "*****   (!SetCamera in FRBeginModeling())" << G4endl;
#endif
		SendStr( FR_SET_CAMERA );

		//----- open device
#if defined DEBUG_FR_SCENE
		G4cerr << "*****   (!OpenDevice in FRBeginModeling())" << G4endl;
#endif
		SendStr( FR_OPEN_DEVICE      );

		//----- begin sending primitives
#if defined DEBUG_FR_SCENE
		G4cerr << "*****   (!BeginModeling in FRBeginModeling())" << G4endl;
#endif
		SendStr( FR_BEGIN_MODELING );  FRflag_in_modeling = true ;

	} // if

} 


/////////////////////////////////////////
// Common to DAWN and DAWNFILE drivers //
/////////////////////////////////////////

#define  G4FRSCENE  G4DAWNFILESceneHandler
#include "G4FRSceneFunc.icc"
#undef   G4FRSCENE 

//////////////////////
// static variables //
//////////////////////

	//----- static variables
G4int G4DAWNFILESceneHandler::fSceneIdCount = 0; 

G4int G4DAWNFILESceneHandler::fSceneCount = 0;   
			// num of existing instances

#endif // G4VIS_BUILD_DAWNFILE_DRIVER
