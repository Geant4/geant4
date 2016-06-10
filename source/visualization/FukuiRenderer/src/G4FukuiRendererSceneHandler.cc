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
// $Id: G4FukuiRendererSceneHandler.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// Satoshi TANAKA, Fri Jun 28 11:34:24 JST 1996
// FukuiRenderer scene.


//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//


#define __G_ANSI_C__

// #define DEBUG_FR_SCENE

     //----- header files
#include <fstream>
#include <string.h>
#include "globals.hh"
#include "G4VisManager.hh"
#include "G4FRConst.hh"
#include "G4FukuiRenderer.hh"
#include "G4FukuiRendererSceneHandler.hh"
#include "G4FukuiRendererViewer.hh"
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
const char  FR_ENV_MULTI_WINDOW[]            = "G4DAWN_MULTI_WINDOW" ;

///////////////////////////
// Driver-dependent part //
///////////////////////////


	//----- G4FukuiRendererSceneHandler, constructor
G4FukuiRendererSceneHandler::G4FukuiRendererSceneHandler (G4FukuiRenderer& system,
					    const G4String& name):
G4VSceneHandler  (system, fSceneIdCount++, name) ,
fSystem   (system)                  ,
fPrimDest (system.GetPrimDest() )   ,
FRflag_in_modeling     (false)      ,
flag_saving_g4_prim    (false)      ,
COMMAND_BUF_SIZE       (G4FRClientServer::SEND_BUFMAX),
fPrec (9), fPrec2 (16)
{

		//----- Connection to FukuiRenderer is set in the first scene
	if( !fSystem.IsConnected() ) 
	{
		if ( getenv( FR_ENV_NAMED_PIPE_CONNECTION ) != NULL &&\
		     strcmp( getenv( FR_ENV_NAMED_PIPE_CONNECTION ), "0" ) )
		{ 
				// Invoke DAWN locally and make connection
				// via named pipe (not supported in AIX etc)
			fSystem.UseBSDUnixDomainAuto();
		} else if( getenv( FR_ENV_SERVER_HOST_NAME ) == NULL  ) 
		{
				// Invoke DAWN locally and make connection
				// via socket
			fSystem.UseInetDomainAuto();
		} else {
				// Connect to remote DAWN via socket
			fSystem.UseInetDomain();
		}
	}

		//----- precision control
	if( getenv( "G4DAWN_PRECISION" ) != NULL ) {
		sscanf( getenv("G4DAWN_PRECISION"), "%d", &fPrec ) ;
	} else {
                fPrec = 9 ;
	}
	fPrec2 = fPrec + 7 ;

} // G4FukuiRendererSceneHandler, constructor


	//----- G4FukuiRendererSceneHandler, destructor
G4FukuiRendererSceneHandler::~G4FukuiRendererSceneHandler () 
{
#if defined DEBUG_FR_SCENE
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** ~G4FukuiRendererSceneHandler" << G4endl;
#endif 
}

//-----
void G4FukuiRendererSceneHandler::FRBeginModeling( void )
{
	if( !FRIsInModeling() )  	
	{
#if defined DEBUG_FR_SCENE
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	    G4cout << "***** G4FukuiRendererSceneHandler::FRBeginModeling (called & started)" << G4endl;
#endif

	  //----- Begin Saving g4.prim file
	  // open g4.prim, ##
	  BeginSavingG4Prim();

	  //----- Send Bounding Box
	  // /BoundingBox
	  SendBoundingBox();

	  //----- drawing device
	  // !Device 
	  // (Note: Not saved in g4.prim) 
	  if( ( getenv( FR_ENV_MULTI_WINDOW ) != NULL      )   && \
	      ( strcmp( getenv( FR_ENV_MULTI_WINDOW ),"0"  )      )  )
	    {
	      ((G4FukuiRendererViewer*)fpViewer)->SendDevice( G4FukuiRendererViewer::FRDEV_XWIN ) ; 
	    } else {
	      ((G4FukuiRendererViewer*)fpViewer)->SendDevice( G4FukuiRendererViewer::FRDEV_PS ) ; 
	    }

	  //----- drawing style
	  // /WireFrame, /Surface etc
	  // (Note: Not saved in g4.prim) 
	  ((G4FukuiRendererViewer*)fpViewer)->SendDrawingStyle() ; 

	  //----- set view parameters
	  //   /CameraPosition, /TargetPoint, 
	  //   /ZoomFactor, /FocalDistance, 
	  //   !GraphicalUserInterface
	  ((G4FukuiRendererViewer*)fpViewer)->SendViewParameters(); 

	  //----- send SET_CAMERA command 
	  //   !SetCamera
#if defined DEBUG_FR_SCENE
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		G4cout << "*****   (!SetCamera in FRBeginModeling())" << G4endl;
#endif
	  SendStr( FR_SET_CAMERA );

	  //----- open device
	  //   !OpenDevice
#if defined DEBUG_FR_SCENE
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		G4cout << "*****   (!OpenDevice in FRBeginModeling())" << G4endl;
#endif
	  SendStr( FR_OPEN_DEVICE      );

	  //----- begin sending primitives
	  //   !BeginModeling
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

#define  G4FRSCENEHANDLER  G4FukuiRendererSceneHandler
#include "G4FRSceneFunc.icc"
#undef   G4FRSCENEHANDLER

//////////////////////
// static variables //
//////////////////////

	//----- static variables
G4int G4FukuiRendererSceneHandler::fSceneIdCount = 0; 

#endif // G4VIS_BUILD_DAWN_DRIVER
