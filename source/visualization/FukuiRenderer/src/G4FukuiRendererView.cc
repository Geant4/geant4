// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FukuiRendererView.cc,v 1.1 1999-01-07 16:14:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Satoshi TANAKA, Fri Jun 28 12:09:11 JST 1996
// FukuiRenderer view - opens window, hard copy, etc.


//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//

#define __G_ANSI_C__
#define G4FukuiRenderer_STRUCTURE_PRIORITY  1.

// #define DEBUG_FR_VIEW

#include "G4ios.hh"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "G4SceneData.hh"
#include "G4Vector3D.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#include "G4FRConst.hh"
#include "G4FukuiRenderer.hh"
#include "G4FukuiRendererScene.hh"
#include "G4FukuiRendererView.hh"


	//----- constants
const char  FR_ENV_MULTI_WINDOW[] = "G4DAWN_MULTI_WINDOW" ;

	//----- G4FukuiRendererView, constructor
G4FukuiRendererView::G4FukuiRendererView (G4FukuiRendererScene& scene,
					  const G4String& name): 
  G4VView (scene, scene.IncrementViewCount (), name), fScene (scene)
{}

	//----- G4FukuiRendererView, destructor
G4FukuiRendererView::~G4FukuiRendererView () 
{}

	//----- G4FukuiRendererView::SetView () 
void G4FukuiRendererView::SetView () 
{
#if defined DEBUG_FR_VIEW
      G4cerr << "***** G4FukuiRendererView::SetView()\n";
#endif 
// Do nothing, since DAWN is running as a different process.
// SendViewParameters () will do this job instead.
}




	//----- G4FukuiRendererView::ClearView()
void
G4FukuiRendererView::ClearView( void )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::ClearView () " << endl;
#endif

		//----- Begin saving data to g4.prim
	fScene.BeginSavingG4Prim();

		//----- Clear old data
	fScene.SendStr( FR_CLEAR_DATA );

}


	//----- G4FukuiRendererView::DrawView () 
void G4FukuiRendererView::DrawView () 
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::DrawView () " << endl;
#endif
		//----- Error recovery
	if( fScene.IsInModeling() )  {
	   G4cerr << "WARNING from FukuiRenderer (DAWN) driver:" << endl;
	   G4cerr << "  You've invoked a drawing command before " << endl;
	   G4cerr << "  completing your previous visualization." << endl;
//	   G4cerr << "  You should have invoked command, /vis~/show/view," << endl;
//	   G4cerr << "  to complete it."   << endl;

	   FlushView();
	}

		//----- Begin Saving g4.prim file
	fScene.BeginSavingG4Prim();

		//----- drawing device
	if( ( getenv( FR_ENV_MULTI_WINDOW ) != NULL      )   && \
	    ( strcmp( getenv( FR_ENV_MULTI_WINDOW ),"0"  )      )  )
	{
		SendDevice( FRDEV_XWIN ) ; 
	} else {
		SendDevice( FRDEV_PS ) ; 
	}

		//----- drawing style
	SendDrawingStyle() ; 

		//----- set view if necessary
	SendViewParameters(); 
		// Camera data should be sent at every time, even if
		// the setting is unchanged at GEANT 4 side.

		//----- Always visit G4 kernel 
	NeedKernelVisit ();
	                           
		//----- Draw
	ProcessView () ;

} // G4FukuiRendererView::DrawView () 



	//----- G4FukuiRendererView::ShowView()
void G4FukuiRendererView::ShowView( void )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::ShowView () " << endl;
#endif

	if( fScene.IsInModeling() ) // if( fScene.flag_in_modeling ) 
	{

			//----- End of Data		
		fScene.FREndModeling();

			//----- Draw all
		fScene.SendStr( FR_DRAW_ALL );

			//----- Close device
		fScene.SendStr( FR_CLOSE_DEVICE );

			//----- End saving data to g4.prim
		fScene.EndSavingG4Prim()              ;

			//----- Wait user clicks drawing Area
		this->Wait();
	}

} // G4FukuiRendererView::ShowView()


	//----- G4FukuiRendererView::FlushView()
void G4FukuiRendererView::FlushView( void )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::Flush () " << endl;
#endif

	if( fScene.IsInModeling() ) // if( fScene.flag_in_modeling ) 
	{

			//----- End of Data		
		fScene.FREndModeling();

			//----- Draw all
		fScene.SendStr( FR_DRAW_ALL );

			//----- Close device
		fScene.SendStr( FR_CLOSE_DEVICE );

			//----- End saving data to g4.prim
		fScene.EndSavingG4Prim()              ;

	}

} // G4FukuiRendererView::FlushView()



	//----- G4FukuiRendererView::Wait()
void
G4FukuiRendererView::Wait()
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::Wait () " << endl;
#endif
  fScene.SendStr    ( FR_WAIT );
  fScene.GetPrimDest().WaitSendBack( FR_WAIT );
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::Wait () : end" << endl;
#endif

}


	//----- G4FukuiRendererView::SendDevice()
void
G4FukuiRendererView::SendDevice( FRDEV dev )
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::SendDevice() " << endl;
#endif


  //	enum {PS=1, XWIN=2, PS2=3, XWIN2=4, OPEN_GL=5, DEVICE_END=6};
  
	if( dev >= FRDEV_PS || dev < FRDEV_DEVICE_END ) {
		fScene.SendStrInt ( FR_DEVICE, dev );
	}
}


	//----- G4FukuiRendererView::SendDrawingStyle()  
void  G4FukuiRendererView::SendDrawingStyle() 
{
#if defined DEBUG_FR_VIEW
	G4cerr << "***** G4FukuiRendererView::SendDrawingStyle() " << endl;
#endif

	G4int  style = fVP.GetDrawingStyle();

	switch( style )
	{
	  case G4ViewParameters::wireframe: 
		fScene.SendStr( FR_WIREFRAME );
		break;
	  case G4ViewParameters::hlr:
		fScene.SendStr( FR_LINES     );
		break;
	  case G4ViewParameters::hsr:
	  case G4ViewParameters::hlhsr:
		fScene.SendStr( FR_SURFACE   );
		break;
	  default:
		fScene.SendStr( FR_WIREFRAME );
		break;
	}

} // G4FukuiRendererView::SendDrawingStyle()  



	//----- G4FukuiRendererView::SendViewParameters () 
void G4FukuiRendererView::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
      G4cerr << "***** G4FukuiRendererView::SendViewParameters()\n";
#endif 

		//----- Magic number to decide camera distance automatically
	const    G4double        HOW_FAR            = 1000.0       ; // to define "infinity"
	const    G4double        MIN_HALF_ANGLE     = 0.01         ;
	const    G4double        MAX_HALF_ANGLE     = 0.499 * M_PI ;

		//----- Send Bounding Box
	fScene.SendBoundingBox();

		//----- (2A) CALC camera distance
		//..... Note: Camera cannot enter inside object
	G4double  camera_distance ;
	G4double  radius = fScene.GetSceneData().GetExtent().GetExtentRadius();

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
		G4cerr << "WARNING from FukuiRenderer (DAWN) driver:" << endl;
		G4cerr << "  Camera cannot enter inside objects"      << endl;
		camera_distance = radius ; 
	}

		//----- (3A) CALC camera direction
	const G4Vector3D& camera_direction \
	  = fVP.GetViewpointDirection().unit();
	const G4double v_angle =  (180.0 / M_PI) * camera_direction.theta() ;
	const G4double h_angle =  (180.0 / M_PI) * camera_direction.phi  () ;

		//----- (2B), (3B) SEND camera position
	fScene.SendStrDouble3( FR_CAMERA_POSITION, 
			       camera_distance, 
			       v_angle, 
			       h_angle ); 

		//----- (4A) CALC target point
	const G4Point3D&  target_point = fVP.GetCurrentTargetPoint();

		//----- (4B) SEND target point
	fScene.SendStrDouble3( FR_TARGET_POINT, 
			       target_point.x(), 
			       target_point.y(), 
			       target_point.z() );

		//----- (5A) CALC zoom factor
	const G4double   zoom_factor  = fVP.GetZoomFactor();

		//----- (5B) SEND zoom factor or focal length
	if( half_view_angle < MIN_HALF_ANGLE ) {

		const G4Point3D&  std_target_point \
	  		= fScene.GetSceneData().GetStandardTargetPoint();

		fScene.SendStrDouble4( FR_ZOOM_FACTOR, 
				       zoom_factor ,
				       std_target_point.x(), 
				       std_target_point.y(), 
				       std_target_point.z());
			// Note that target point, camera position, 
			// and bounding box have already been sent above.
			// The std_target_point is necessary to
			// Calc focal distance from the zoom factor.
	} else {
		const G4double FR_HALF_SCREEN_SIZE = 0.5 ;
		G4double  focal_distance \
		  = FR_HALF_SCREEN_SIZE / tan( half_view_angle ); 
		focal_distance *= zoom_factor ;
		fScene.SendStrDouble ( FR_FOCAL_DISTANCE, focal_distance );
	}

		//----- INVOKE GUI: not executed in the default setting
	if( fScene.GetSystem().IsGUIMode() ) {
			//----- send GUI command
		fScene.SendStr( FR_GUI );

			//----- wait the same command is sent back:
			//..... This avoids to send many data before
			//..... GUI session is over.
		fScene.GetPrimDest().WaitSendBack( FR_GUI );
	}

		//----- SET CAMERA
	fScene.SendStr( FR_SET_CAMERA );

} // G4FukuiRendererView::SendViewParameters () 

#endif // G4VIS_BUILD_DAWN_DRIVER
