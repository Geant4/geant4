// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DAWNFILEViewer.cc,v 1.2 1999-01-11 00:47:21 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA
// DAWNFILE view - opens window, hard copy, etc.


//=================//
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
//=================//

#define __G_ANSI_C__
#define G4DAWNFILE_STRUCTURE_PRIORITY  1.

// #define DEBUG_FR_VIEW

#include "G4ios.hh"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "G4Scene.hh"
#include "G4Vector3D.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#include "G4FRConst.hh"
#include "G4DAWNFILE.hh"
#include "G4DAWNFILESceneHandler.hh"
#include "G4DAWNFILEViewer.hh"



	//----- constants
const char  FR_ENV_MULTI_WINDOW[] = "G4DAWN_MULTI_WINDOW" ;

	//----- G4DAWNFILEViewer, constructor
G4DAWNFILEViewer::G4DAWNFILEViewer (G4DAWNFILESceneHandler& scene,
				const G4String& name): 
  G4VViewer (scene, scene.IncrementViewCount (), name), fSceneHandler (scene)
{
	// Set a g4.prim-file viewer 
	strcpy( fG4PrimViewer, "dawn" ); 
	if( getenv( "G4DAWNFILE_VIEWER" ) != NULL ) {
		strcpy( fG4PrimViewer, getenv( "G4DAWNFILE_VIEWER" ) ) ;			
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
	strcpy( fPSViewer, "ghostview" ); 
	if( getenv( "G4DAWNFILE_PS_VIEWER" ) != NULL ) {
		strcpy( fPSViewer, getenv( "G4DAWNFILE_PS_VIEWER" ) ) ;			
	} 

}

	//----- G4DAWNFILEViewer, destructor
G4DAWNFILEViewer::~G4DAWNFILEViewer () 
{}

	//----- G4DAWNFILEViewer::SetView () 
void G4DAWNFILEViewer::SetView () 
{
#if defined DEBUG_FR_VIEW
      G4cerr << "***** G4DAWNFILEViewer::SetView()\n";
#endif 
// Do nothing, since DAWN is running as a different process.
// SendViewParameters () will do this job instead.
}


	//----- G4DAWNFILEViewer::ClearView()
void
G4DAWNFILEViewer::ClearView( void )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4DAWNFILEViewer::ClearView () " << endl;
#endif

		//----- Begin saving data to g4.prim
	fSceneHandler.BeginSavingG4Prim();

		//----- Clear old data
	fSceneHandler.SendStr( FR_CLEAR_DATA );

}


	//----- G4DAWNFILEViewer::DrawView () 
void G4DAWNFILEViewer::DrawView () 
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4DAWNFILEViewer::DrawView () " << endl;
#endif
		//----- Error recovery
	if( fSceneHandler.IsInModeling() )  {
	   G4cerr << "WARNING from DAWNFILE driver:" << endl;
	   G4cerr << "  You've invoked a drawing command before " << endl;
	   G4cerr << "  completing your previous visualization." << endl;
//	   G4cerr << "  You should have invoked command, /vis~/show/view," << endl;
//	   G4cerr << "  to complete it."   << endl;

	   FlushView();
	}

		//----- Begin Saving g4.prim file
	fSceneHandler.BeginSavingG4Prim();

		//----- set view if necessary
	SendViewParameters(); 
		// Camera data should be sent at every time, even if
		// the setting is unchanged at GEANT 4 side.

		//----- Always visit G4 kernel 
	NeedKernelVisit ();
	                           
		//----- Draw
	ProcessView () ;

} // G4DAWNFILEViewer::DrawView () 



	//----- G4DAWNFILEViewer::ShowView()
void G4DAWNFILEViewer::ShowView( void )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4DAWNFILEViewer::ShowView () " << endl;
#endif

	if( fSceneHandler.IsInModeling() ) // if( fSceneHandler.flag_in_modeling ) 
	{

			//----- End of Data		
		fSceneHandler.FREndModeling();

			//----- Draw all
		fSceneHandler.SendStr( FR_DRAW_ALL );

			//----- Close device
		fSceneHandler.SendStr( FR_CLOSE_DEVICE );

			//----- End saving data to g4.prim
		fSceneHandler.EndSavingG4Prim()              ;


			//----- string for viewer invocation
		if ( !strcmp( fG4PrimViewer, "NONE" ) ) {
		
			strcpy( fG4PrimViewerInvocation, "" );
		} else {

			strcpy( fG4PrimViewerInvocation, fG4PrimViewer );
			strcat( fG4PrimViewerInvocation, " ");
			strcat( fG4PrimViewerInvocation, fSceneHandler.GetG4PrimFileName() );
		}


		//----- Invoke DAWN
		if( false == G4FRofstream::DoesFileExist( fSceneHandler.GetG4PrimFileName() ) )   
		{
			G4cout << "ERROR: Failed to generate file  ";
			G4cout << fSceneHandler.GetG4PrimFileName() << endl;

		} else 	if( strcmp( GetG4PrimViewerInvocation(), "" ) )  
		{
			G4cout << "File  " << fSceneHandler.GetG4PrimFileName() ;
			G4cout << "  is generated." << endl;
			G4cout << GetG4PrimViewerInvocation() << endl;
			system( GetG4PrimViewerInvocation() );

		} else { // no view, i.e., only file generation
			G4cout << "File  " << fSceneHandler.GetG4PrimFileName() ; 
			G4cout << "  is generated." << endl;
			G4cout << "No viewer is invoked." << endl;
		}

	}

} // G4DAWNFILEViewer::ShowView()


	//----- G4DAWNFILEViewer::FlushView()
void G4DAWNFILEViewer::FlushView( void )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4DAWNFILEViewer::Flush () " << endl;
#endif

	ShowView();
} 



	//----- G4DAWNFILEViewer::SendDrawingStyleToDAWNGUI( ostream& out ) 
void  G4DAWNFILEViewer::SendDrawingStyleToDAWNGUI( ostream& out ) 
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4DAWNFILEViewer::SendDrawingStyle() " << endl;
#endif

	G4int  style = fVP.GetDrawingStyle();

	enum {	FR_WIREFRAME = 1, FR_WF_STORED = 2, FR_HID =3  , \
		FR_HID2      = 4, FR_HID3      = 5, FR_DRAWING_MODE_END = 6 };

	switch( style )
	{
	  case G4ViewParameters::wireframe: 
		out <<  FR_WIREFRAME << endl;
		break;
	  case G4ViewParameters::hlr:
		out <<  FR_HID2      << endl; // LINE
		break;
	  case G4ViewParameters::hsr:
	  case G4ViewParameters::hlhsr:
		out <<  FR_HID       << endl; // SURFACE
		break;
	  default:
		out <<  FR_WIREFRAME << endl;
		break;
	}

} // G4DAWNFILEViewer::SendDrawingStyle()  



	//----- G4DAWNFILEViewer::SendViewParameters () 
void G4DAWNFILEViewer::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
      G4cerr << "***** G4DAWNFILEViewer::SendViewParameters()\n";
#endif 

		//----- Magic number to decide camera distance automatically
	const    G4double        HOW_FAR            = 1000.0       ; // to define "infinity"
	const    G4double        MIN_HALF_ANGLE     = 0.01         ;
	const    G4double        MAX_HALF_ANGLE     = 0.499 * M_PI ;

		//----- Send Bounding Box
	fSceneHandler.SendBoundingBox();

		//----- CALC camera distance
		//..... Note: Camera cannot enter inside object
	G4double  camera_distance ;
	G4double  radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();

	G4double half_view_angle  = fabs ( fVP.GetFieldHalfAngle () ) ;
	if( half_view_angle > MAX_HALF_ANGLE ) { 
	  half_view_angle = MAX_HALF_ANGLE ; 
	} 

	if( half_view_angle < MIN_HALF_ANGLE ) {
			//----- infinity (or ortho projection)
		camera_distance = radius * HOW_FAR ;  
	} else {
			//----- Calc camera distance from half view angle
		camera_distance = radius / sin ( half_view_angle );
		camera_distance -= fVP.GetDolly();
	}

	if ( camera_distance < radius ) { 
		G4cerr << "WARNING from DAWNFILE driver:" << endl;
		G4cerr << "  Camera cannot enter inside objects"      << endl;
		camera_distance = radius ; 
	}

		//----- CALC camera direction
	const G4Vector3D& camera_direction \
	  = fVP.GetViewpointDirection().unit();
	const G4double v_angle =  (180.0 / M_PI) * camera_direction.theta() ;
	const G4double h_angle =  (180.0 / M_PI) * camera_direction.phi  () ;


	//########### Generation of the file .DAWN.history for DAWN GUI
		//-----
	ofstream gui_out (".DAWN.history") ; 

	// ######### P1 

		//----- camera position
	gui_out << camera_distance << endl;
	gui_out << v_angle  << endl ;
        gui_out << h_angle  << endl ;
        gui_out << "0"  << endl     ; // auto target

		//----- target point 
	const G4Point3D&  target_point = fVP.GetCurrentTargetPoint();
	gui_out << target_point.x()          << endl ;
	gui_out << target_point.y()          << endl ;
	gui_out << target_point.z()          << endl ;

		//----- Magnification
	const G4double   zoom_factor  = fVP.GetZoomFactor();
	if( half_view_angle < MIN_HALF_ANGLE ) {

		gui_out << zoom_factor << endl;

	} else {
		const G4double FR_HALF_SCREEN_SIZE = 0.5 ;
		G4double  focal_distance \
		  = FR_HALF_SCREEN_SIZE / tan( half_view_angle ); 
		focal_distance *= zoom_factor ;

		gui_out << "fd" << focal_distance << endl;

	}
	SendDrawingStyleToDAWNGUI( gui_out ) ; // gui_out, viewing mode
	gui_out << "0.001" << endl           ; // 3D Tolerance 
	gui_out << "0"     << endl           ; // not display parameters


	// ######### P2
	gui_out << 1 << endl;   // Source light 
	gui_out << 1 << endl;   
	gui_out << 1 << endl;   
	gui_out << 0.5 << endl; // Ambient light
	gui_out << 0.5 << endl;
	gui_out << 0.5 << endl;
	gui_out << 19.0 << endl; // Light direction (Polar)
	gui_out << 71.0 << endl; // Light direction (Azimuthal)

	// ######### P3
	gui_out << 0.1 << endl;    // Real edge width
	gui_out << 0.1 << endl;    // outline   width
	gui_out << 0.1 << endl;    // aux edge  width
	gui_out << 3   << endl;      // aux edge  style
	gui_out << 70.0<< endl;   // aux-edge threshold angle
	gui_out << 0.1 << endl;       // line width
	gui_out << 0   << endl;        // haloing
	gui_out << 1   << endl;        // Dashed edged for back faces

	//######### P4
               //----- drawing device
  //	enum {PS=1, XWIN=2, PS2=3, XWIN2=4, OPEN_GL=5, DEVICE_END=6};
        if( ( getenv( FR_ENV_MULTI_WINDOW ) != NULL      )   && \
            ( strcmp( getenv( FR_ENV_MULTI_WINDOW ),"0"  )      )  )
        {
                gui_out << 2 << endl; // OpenWindow
        } else {
                gui_out << 1 << endl; // Invoke PS viewer
        }
	gui_out << GetPSViewer() << endl; // PS viewer
	gui_out << 0;         // Add showpage 

	

	gui_out.close();
	//########### end of generating file .DAWN.history 


		//----- SET CAMERA
	fSceneHandler.SendStr( FR_SET_CAMERA );


} // G4DAWNFILEViewer::SendViewParameters () 


#endif // G4VIS_BUILD_DAWNFILE_DRIVER
