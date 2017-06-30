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
// $Id: G4DAWNFILESceneHandler.cc 104015 2017-05-08 07:28:08Z gcosmo $
//
// Satoshi TANAKA
// DAWNFILE scene.


#define __G_ANSI_C__

// #define DEBUG_FR_SCENE

     //----- header files
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <iomanip>
#include "globals.hh"
#include "G4VisManager.hh"
#include "G4FRConst.hh"
#include "G4DAWNFILE.hh"
#include "G4DAWNFILESceneHandler.hh"
#include "G4DAWNFILEViewer.hh"
#include "G4Point3D.hh"
#include "G4VisAttributes.hh"
#include "G4Scene.hh"
#include "G4Transform3D.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyline.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4VPhysicalVolume.hh"

//----- constants
const char  FR_ENV_CULL_INVISIBLE_OBJECTS [] = "G4DAWN_CULL_INVISIBLE_OBJECTS";
const char  G4PRIM_FILE_HEADER      [] = "g4_";
const char  DEFAULT_G4PRIM_FILE_NAME[] = "g4_0000.prim";

// const int   FR_MAX_FILE_NUM = 1 ;
// const int   FR_MAX_FILE_NUM = 5 ;
// const int   FR_MAX_FILE_NUM = 10 ;
// const int   FR_MAX_FILE_NUM = 15 ;
// const int   FR_MAX_FILE_NUM = 20 ;
   const int   FR_MAX_FILE_NUM = 100 ;


///////////////////////////
// Driver-dependent part //
///////////////////////////


	//----- G4DAWNFILESceneHandler, constructor
G4DAWNFILESceneHandler::G4DAWNFILESceneHandler (G4DAWNFILE& system, const G4String& name):
G4VSceneHandler  (system, fSceneIdCount++, name) ,
fSystem   (system)                        ,
fPrimDest ()                              ,
FRflag_in_modeling     (false)            ,
flag_saving_g4_prim    (false)            ,
COMMAND_BUF_SIZE       (G4FRofstream::SEND_BUFMAX),
fPrec (9), fPrec2 (16)
{
	// g4.prim filename and its directory
	if ( getenv( "G4DAWNFILE_DEST_DIR" ) == NULL ) {
		strcpy( fG4PrimDestDir , "" )                      ;  // output dir
		strcpy( fG4PrimFileName, DEFAULT_G4PRIM_FILE_NAME );  // filename
	} else {
		strcpy( fG4PrimDestDir , getenv( "G4DAWNFILE_DEST_DIR" ) ); // output dir
		strcpy( fG4PrimFileName, DEFAULT_G4PRIM_FILE_NAME        ); // filename 
	}
		
	// maximum number of g4.prim files in the dest directory
	fMaxFileNum = FR_MAX_FILE_NUM ; // initialization
	if ( getenv( "G4DAWNFILE_MAX_FILE_NUM" ) != NULL ) {	
		
		sscanf( getenv("G4DAWNFILE_MAX_FILE_NUM"), "%d", &fMaxFileNum ) ;

	} else {
		fMaxFileNum = FR_MAX_FILE_NUM ;
	}
	if( fMaxFileNum < 1 ) { fMaxFileNum = 1 ; }


		//----- precision control
	if( getenv( "G4DAWNFILE_PRECISION" ) != NULL ) {
		sscanf( getenv("G4DAWNFILE_PRECISION"), "%d", &fPrec ) ;
	} else {
                fPrec = 9 ;
	}
	fPrec2 = fPrec + 7 ;

} 


	//----- G4DAWNFILESceneHandler, destructor
G4DAWNFILESceneHandler::~G4DAWNFILESceneHandler () 
{
#if defined DEBUG_FR_SCENE
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** ~G4DAWNFILESceneHandler" << G4endl;
#endif 
	if (fPrimDest.IsOpen()) 
	  {
			//----- End of modeling
			// !EndModeling, !DrawAll, !CloseDevice,
			// close g4.prim
		FREndModeling();
	}
}

//-----
void	G4DAWNFILESceneHandler::SetG4PrimFileName() 
{
	// g4_0000.prim, g4_0001.prim, ..., g4_MAX_FILE_INDEX.prim
	const int MAX_FILE_INDEX = fMaxFileNum - 1 ;

	// dest directory (null if no environmental variables is set)
	strcpy ( fG4PrimFileName, fG4PrimDestDir) ; 

	// create full path name (default)
	strcat ( fG4PrimFileName, DEFAULT_G4PRIM_FILE_NAME );

	// Automatic updation of file names
	for( int i = 0 ; i < fMaxFileNum ; i++) { 

		// Message in the final execution
		if( i == MAX_FILE_INDEX ) 
		{
		  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
		    G4cout << "==========================================="   << G4endl; 
		    G4cout << "WARNING MESSAGE from DAWNFILE driver:      "   << G4endl;
		    G4cout << "  This file name is the final one in the   "   << G4endl;
		    G4cout << "  automatic updation of the output file name." << G4endl; 
		    G4cout << "  You may overwrite existing files, i.e.   "   << G4endl; 
		    G4cout << "  g4_XXXX.prim and g4_XXXX.eps             "   << G4endl;
		    G4cout << "==========================================="   << G4endl; 
		  }
		}

		// re-determine file name as G4DAWNFILE_DEST_DIR/g4_XXXX.prim
		std::ostringstream filename; filename
		<< fG4PrimDestDir << G4PRIM_FILE_HEADER
		<< std::setw(4) << std::setfill('0') << i << ".prim";
		strncpy(fG4PrimFileName,filename.str().c_str(),sizeof(fG4PrimFileName));

		// check validity of the file name
		std::ifstream  fin ; 
		fin.open(fG4PrimFileName) ;
		if(!fin) { 
			// new file	
			fin.close();  
			break; 
		} else { 
			// already exists (try next) 
			fin.close(); 
		} 

	} // for 

	G4cout << "===========================================    " << G4endl; 
	G4cout << "Output file: " <<    fG4PrimFileName             << G4endl; 
	G4cout << "Destination directory (current dir if NULL): "       << fG4PrimDestDir    << G4endl; 
	G4cout << "Maximal number of files in the destination directory: " << fMaxFileNum << G4endl; 
	G4cout << "Note:                                                " << G4endl; 
	G4cout << "  * The maximal number is customizable as:           " << G4endl;
	G4cout << "       % setenv  G4DAWNFILE_MAX_FILE_NUM  number " << G4endl;        
	G4cout << "  * The destination directory is customizable as:" << G4endl;
	G4cout << "       % setenv  G4DAWNFILE_DEST_DIR  dir_name/  " << G4endl;        
	G4cout << "     ** Do not forget \"/\" at the end of the    " << G4endl;              
	G4cout << "        dir_name, e.g. \"./tmp/\".  " << G4endl;              
	G4cout << "===========================================      " << G4endl; 

} // G4DAWNFILESceneHandler::SetG4PrimFileName()


//-----
void	G4DAWNFILESceneHandler::BeginSavingG4Prim( void ) 
{
#if defined DEBUG_FR_SCENE
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** BeginSavingG4Prim (called)\n";
#endif

	if( !IsSavingG4Prim() ) 
	{ 
#if defined DEBUG_FR_SCENE
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	        G4cout << "*****                   (started) " ;
	        G4cout << "(open g4.prim, ##)"  << G4endl;
	  }
#endif
		SetG4PrimFileName() ; // result set to fG4PrimFileName
		fPrimDest.Open(fG4PrimFileName)   ;

		SendStr( FR_G4_PRIM_HEADER   )    ; 
		flag_saving_g4_prim = true        ; 
	} 
}

void	G4DAWNFILESceneHandler::EndSavingG4Prim  ( void ) 
{
#if defined DEBUG_FR_SCENE
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** EndSavingG4Prim (called)\n";
#endif

	if(  IsSavingG4Prim() )
	{
#if defined DEBUG_FR_SCENE
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	    G4cout << "*****                 (started) (close g4.prim)" << G4endl;
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
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	    G4cout << "***** G4DAWNFILESceneHandler::FRBeginModeling (called & started)" << G4endl;
#endif

			//----- Send saving command and heading comment
		BeginSavingG4Prim();

			//----- Send bounding box command
		SendBoundingBox();

			//----- send SET_CAMERA command 
#if defined DEBUG_FR_SCENE
		if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		  G4cout << "*****   (!SetCamera in FRBeginModeling())" << G4endl;
#endif
		SendStr( FR_SET_CAMERA );

		//----- open device
#if defined DEBUG_FR_SCENE
		if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		  G4cout << "*****   (!OpenDevice in FRBeginModeling())" << G4endl;
#endif
		SendStr( FR_OPEN_DEVICE      );

		//----- begin sending primitives
#if defined DEBUG_FR_SCENE
		if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		  G4cout << "*****   (!BeginModeling in FRBeginModeling())" << G4endl;
#endif
		SendStr( FR_BEGIN_MODELING );  FRflag_in_modeling = true ;

	} // if

} 


/////////////////////////////////////////
// Common to DAWN and DAWNFILE drivers //
/////////////////////////////////////////

#define  G4FRSCENEHANDLER  G4DAWNFILESceneHandler
#include "G4FRSceneFunc.icc"
#undef   G4FRSCENEHANDLER 

//////////////////////
// static variables //
//////////////////////

	//----- static variables
G4int G4DAWNFILESceneHandler::fSceneIdCount = 0; 
