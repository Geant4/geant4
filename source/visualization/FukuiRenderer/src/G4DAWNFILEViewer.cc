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
//
// Satoshi TANAKA
// DAWNFILE view - opens window, hard copy, etc.


#define __G_ANSI_C__
#define G4DAWNFILE_STRUCTURE_PRIORITY  1.

// #define DEBUG_FR_VIEW

#include "G4ios.hh"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "G4VisManager.hh"
#include "G4Scene.hh"
#include "G4Vector3D.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PhysicalConstants.hh"

#include "G4FRConst.hh"
#include "G4DAWNFILE.hh"
#include "G4DAWNFILESceneHandler.hh"
#include "G4DAWNFILEViewer.hh"



	//----- constants
const char  FR_ENV_MULTI_WINDOW [] = "G4DAWN_MULTI_WINDOW" ;
const char  FR_ENV_MULTI_WINDOW2[] = "G4DAWNFILE_MULTI_WINDOW" ;

	//----- G4DAWNFILEViewer, constructor
G4DAWNFILEViewer::G4DAWNFILEViewer (G4DAWNFILESceneHandler& sceneHandler,
				const G4String& name): 
  G4VViewer (sceneHandler,
	     sceneHandler.IncrementViewCount (),
	     name),
  fSceneHandler (sceneHandler)
{
	// Set a g4.prim-file viewer 
	strcpy( fG4PrimViewer, "dawn" ); 
	if( std::getenv( "G4DAWNFILE_VIEWER" ) != NULL ) {
		strcpy( fG4PrimViewer, std::getenv( "G4DAWNFILE_VIEWER" ) ) ;			
	} 

	// string for viewer invocation
	if ( !strcmp( fG4PrimViewer, "NONE" ) ) {
		
		strcpy( fG4PrimViewerInvocation, "" );
	} else {

		strcpy( fG4PrimViewerInvocation, fG4PrimViewer );
		strcat( fG4PrimViewerInvocation, " ");
		strcat( fG4PrimViewerInvocation, fSceneHandler.GetG4PrimFileName() );
	}

	// Set a PostScript Viewer
//	strcpy( fPSViewer, "ghostview" ); 
	strcpy( fPSViewer, "gv" ); 
	if( std::getenv( "G4DAWNFILE_PS_VIEWER" ) != NULL ) {
		strcpy( fPSViewer, std::getenv( "G4DAWNFILE_PS_VIEWER" ) ) ;			
	} 

}

	//----- G4DAWNFILEViewer, destructor
G4DAWNFILEViewer::~G4DAWNFILEViewer () 
{}

	//----- G4DAWNFILEViewer::SetView () 
void G4DAWNFILEViewer::SetView () 
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
    G4cout << "***** G4DAWNFILEViewer::SetView(): No effects" << G4endl;
#endif 
// Do nothing, since DAWN is running as a different process.
// SendViewParameters () will do this job instead.
}


	//----- G4DAWNFILEViewer::ClearView()
void
G4DAWNFILEViewer::ClearView( void )
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4DAWNFILEViewer::ClearView (): No effects " << G4endl;
#endif
	if (fSceneHandler.fPrimDest.IsOpen()) {
	  fSceneHandler.fPrimDest.Close();
	  // Re-open with same filename...
	  fSceneHandler.fPrimDest.Open(fSceneHandler.fG4PrimFileName);
	  fSceneHandler.SendStr( FR_G4_PRIM_HEADER );
	  fSceneHandler.FRflag_in_modeling = false;
	  fSceneHandler.FRBeginModeling();
	}
}


	//----- G4DAWNFILEViewer::DrawView () 
void G4DAWNFILEViewer::DrawView () 
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4DAWNFILEViewer::DrawView () " << G4endl;
#endif
		//----- 
	fSceneHandler.FRBeginModeling() ;

		//----- Always visit G4 kernel 
	NeedKernelVisit ();
	                           
		//----- Draw
	ProcessView () ;

} // G4DAWNFILEViewer::DrawView () 



	//----- G4DAWNFILEViewer::ShowView()
void G4DAWNFILEViewer::ShowView( void )
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4DAWNFILEViewer::ShowView () " << G4endl;
#endif

	if( fSceneHandler.FRIsInModeling() ) 
	{
			//----- End of modeling
			// !EndModeling, !DrawAll, !CloseDevice,
			// close g4.prim
		fSceneHandler.FREndModeling();

			//----- Output DAWN GUI file 
		SendViewParameters(); 

			//----- string for viewer invocation
		if ( !strcmp( fG4PrimViewer, "NONE" ) ) {
		
			strcpy( fG4PrimViewerInvocation, "" );
		} else {

			strcpy( fG4PrimViewerInvocation, fG4PrimViewer );
			strcat( fG4PrimViewerInvocation, " ");
			strcat( fG4PrimViewerInvocation, fSceneHandler.GetG4PrimFileName() );
		}


		//----- Invoke DAWN
		G4cout << G4endl ;
		if( false == G4FRofstream::DoesFileExist( fSceneHandler.GetG4PrimFileName() ) )   
		{
			G4cout << "ERROR: Failed to generate file  ";
			G4cout << fSceneHandler.GetG4PrimFileName() << G4endl;

		} else 	if( strcmp( GetG4PrimViewerInvocation(), "" ) )  
		{
			G4cout << "File  " << fSceneHandler.GetG4PrimFileName() ;
			G4cout << "  is generated." << G4endl;
			G4cout << GetG4PrimViewerInvocation() << G4endl;
			(void) system( GetG4PrimViewerInvocation() );

		} else { // no view, i.e., only file generation
			G4cout << "File  " << fSceneHandler.GetG4PrimFileName() ; 
			G4cout << "  is generated." << G4endl;
			G4cout << "No viewer is invoked." << G4endl;
		}

	}

} // G4DAWNFILEViewer::ShowView()


	//----- G4DAWNFILEViewer::SendDrawingStyleToDAWNGUI( std::ostream& out ) 
void  G4DAWNFILEViewer::SendDrawingStyleToDAWNGUI( std::ostream& out ) 
{
///////////////////////
//#if defined DEBUG_FR_VIEW
//  G4cout << "***** G4DAWNFILEViewer::SendDrawingStyleToDAWNGUI()" << G4endl;
//#endif
//////////////////////

	G4int  style = fVP.GetDrawingStyle();

	enum {	FR_WIREFRAME_STYLE = 1, FR_WF_STORED_STYLE = 2,
		FR_HID_STYLE       = 3,	FR_HID2_STYLE      = 4,
		FR_HID3_STYLE      = 5, FR_DRAWING_MODE_END_STYLE = 6 };

	switch( style )
	{
	  case G4ViewParameters::wireframe: 
		out <<  FR_WIREFRAME_STYLE << G4endl;
		break;
	  case G4ViewParameters::hlr:
		out <<  FR_HID2_STYLE      << G4endl; // LINE
		break;
	  case G4ViewParameters::hsr:
	  case G4ViewParameters::hlhsr:
		out <<  FR_HID_STYLE       << G4endl; // SURFACE
		break;
	  default:
		out <<  FR_WIREFRAME_STYLE << G4endl;
		break;
	}

} // G4DAWNFILEViewer::SendDrawingStyle()  



//----- 
void G4DAWNFILEViewer::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
    G4cout << "***** G4DAWNFILEViewer::SendViewParameters()  "; 
    G4cout << "(GUI parameters)" << G4endl;
  }
#endif 

		//----- Magic number to decide camera distance automatically
	const    G4double        HOW_FAR            = 1000.0       ; // to define "infinity"
	const    G4double        MIN_HALF_ANGLE     = 0.01         ;
	const    G4double        MAX_HALF_ANGLE     = 0.499 * pi ;

		//----- CALC camera distance
		//..... Note: Camera cannot enter inside object
	G4double  camera_distance ;
	G4double  radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();

	G4double half_view_angle  = std::fabs ( fVP.GetFieldHalfAngle () ) ;
	if( half_view_angle > MAX_HALF_ANGLE ) { 
	  half_view_angle = MAX_HALF_ANGLE ; 
	} 

	if( half_view_angle < MIN_HALF_ANGLE ) {
			//----- infinity (or ortho projection)
		camera_distance = radius * HOW_FAR ;  
	} else {
			//----- Calc camera distance from half view angle
		camera_distance = radius / std::sin ( half_view_angle );
		camera_distance -= fVP.GetDolly();
	}

	if ( camera_distance < radius ) { 
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
		G4cout << "WARNING from DAWNFILE driver:" << G4endl;
		G4cout << "  Camera cannot enter inside objects"      << G4endl;
	  }
		camera_distance = radius ; 
	}

		//----- CALC camera direction
	const G4Vector3D& camera_direction \
	  = fVP.GetViewpointDirection().unit();
	const G4double v_angle =  (180.0 / pi) * camera_direction.theta() ;
	const G4double h_angle =  (180.0 / pi) * camera_direction.phi  () ;


	//########### Generation of the file .DAWN.history for DAWN GUI
		//-----
	std::ofstream gui_out (".DAWN_1.history") ; 

	// ######### P1 

		//----- camera position
	gui_out << camera_distance << G4endl;
	gui_out << v_angle  << G4endl ;
        gui_out << h_angle  << G4endl ;
        gui_out << "0"  << G4endl     ; // auto target

		//----- target point 
	const G4Point3D&  target_point
	  = fSceneHandler.GetScene()->GetStandardTargetPoint()
	  + fVP.GetCurrentTargetPoint();
	gui_out << target_point.x()          << G4endl ;
	gui_out << target_point.y()          << G4endl ;
	gui_out << target_point.z()          << G4endl ;

		//----- Magnification
	const G4double   zoom_factor  = fVP.GetZoomFactor();
	if( half_view_angle < MIN_HALF_ANGLE ) {

		gui_out << zoom_factor << G4endl;

	} else {
		const G4double FR_HALF_SCREEN_SIZE = 0.5 ;
		G4double  focal_distance \
		  = FR_HALF_SCREEN_SIZE / std::tan( half_view_angle ); 
		focal_distance *= zoom_factor ;

		gui_out << "fd" << focal_distance << G4endl;

	}
	SendDrawingStyleToDAWNGUI( gui_out ) ; // gui_out, viewing mode
	gui_out << "0.001" << G4endl           ; // 3D Tolerance 
	gui_out << "0"     << G4endl           ; // not display parameters


	// ######### P2
	gui_out << 1 << G4endl;   // Source light 
	gui_out << 1 << G4endl;   
	gui_out << 1 << G4endl;   
	gui_out << 0.5 << G4endl; // Ambient light
	gui_out << 0.5 << G4endl;
	gui_out << 0.5 << G4endl;
	gui_out << 19.0 << G4endl; // Light direction (Polar)
	gui_out << 71.0 << G4endl; // Light direction (Azimuthal)

	// ######### P3
	gui_out << 0.1 << G4endl;    // Real edge width
	gui_out << 0.1 << G4endl;    // outline   width
	gui_out << 0.1 << G4endl;    // aux edge  width
	gui_out << 3   << G4endl;      // aux edge  style
	gui_out << 70.0<< G4endl;   // aux-edge threshold angle
	gui_out << 0.1 << G4endl;       // line width
	gui_out << 0   << G4endl;        // haloing
	gui_out << 1   << G4endl;        // Dashed edged for back faces

	//######### P4
               //----- drawing device
  //	enum {PS=1, XWIN=2, PS2=3, XWIN2=4, OPEN_GL=5, DEVICE_END=6};
        if( ( ( std::getenv( FR_ENV_MULTI_WINDOW ) != NULL        ) && \
              ( strcmp( std::getenv( FR_ENV_MULTI_WINDOW ),"0"  ) )       ) || \
            ( ( std::getenv( FR_ENV_MULTI_WINDOW2 ) != NULL        ) && \
              ( strcmp( std::getenv( FR_ENV_MULTI_WINDOW2 ),"0"  ) )      )     )
        {
                gui_out << 2 << G4endl; // OpenWindow
        } else {
                gui_out << 1 << G4endl; // Invoke PS viewer
        }

	gui_out << GetPSViewer() << G4endl; // PS viewer
	gui_out << 1 << G4endl            ; // Do not add showpage 
	gui_out << 0 << G4endl            ; // Non-append mode

	gui_out.close();
	//########### end of generating file .DAWN.history 


} 
