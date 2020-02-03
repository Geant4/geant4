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

#include "G4VisManager.hh"
#include "G4Scene.hh"
#include "G4Vector3D.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PhysicalConstants.hh"

#include "G4FRConst.hh"
#include "G4FukuiRenderer.hh"
#include "G4FukuiRendererSceneHandler.hh"
#include "G4FukuiRendererViewer.hh"


//----- Constructor
G4FukuiRendererViewer::G4FukuiRendererViewer (G4FukuiRendererSceneHandler& sceneHandler,
					  const G4String& name): 
  G4VViewer (sceneHandler,
	     sceneHandler.IncrementViewCount (),
	     name),
  fSceneHandler (sceneHandler)
{}

//----- Destructor
G4FukuiRendererViewer::~G4FukuiRendererViewer () 
{}

//----- 
void G4FukuiRendererViewer::SetView () 
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
    G4cout << "***** G4FukuiRendererViewer::SetView(): No effects" << G4endl;
#endif 
// Do nothing, since DAWN is running as a different process.
// SendViewParameters () will do this job instead.
}

//----- 
void
G4FukuiRendererViewer::ClearView( void )
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
    G4cout << "***** G4FukuiRendererViewer::ClearView (): No effects " << G4endl;
#endif

}


//----- 
void G4FukuiRendererViewer::DrawView () 
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4FukuiRendererViewer::DrawView () " << G4endl;
#endif

	//----- Begin modeling 3D 
	fSceneHandler.FRBeginModeling();	

	//----- Always visit G4 kernel 
	NeedKernelVisit ();
	                           
	//----- Draw
	ProcessView () ;

} 


//----- 
void G4FukuiRendererViewer::ShowView( void )
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4FukuiRendererViewer::ShowView () " << G4endl;
#endif

	if( fSceneHandler.FRIsInModeling() ) 
	{
			//----- End of modeling
			// !EndModeling, !DrawAll, !CloseDevice,
			// close g4.prim
		fSceneHandler.FREndModeling();

			//----- Wait user clicks drawing Area
		this->Wait();
	}

} 


//----- 
void  G4FukuiRendererViewer::Wait()
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4FukuiRendererViewer::Wait () : Begin" << G4endl;
#endif
  fSceneHandler.SendStr    ( FR_WAIT );
  fSceneHandler.GetPrimDest().WaitSendBack( FR_WAIT );
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4FukuiRendererViewer::Wait () : end" << G4endl;
#endif

}


//----- 
void
G4FukuiRendererViewer::SendDevice( FRDEV dev )
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4FukuiRendererViewer::SendDevice() " << G4endl;
#endif

  //	enum {PS=1, XWIN=2, PS2=3, XWIN2=4, OPEN_GL=5, DEVICE_END=6};
  
	if( dev >= FRDEV_PS || dev < FRDEV_DEVICE_END ) {
		fSceneHandler.SendStrInt ( FR_DEVICE, dev );
	}
}


//----- 
void  G4FukuiRendererViewer::SendDrawingStyle() 
{
#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** G4FukuiRendererViewer::SendDrawingStyle() " << G4endl;
#endif

	G4int  style = fVP.GetDrawingStyle();

	switch( style )
	{
	  case G4ViewParameters::wireframe: 
		fSceneHandler.SendStr( FR_WIREFRAME );
		break;
	  case G4ViewParameters::hlr:
		fSceneHandler.SendStr( FR_LINES     );
		break;
	  case G4ViewParameters::hsr:
	  case G4ViewParameters::hlhsr:
		fSceneHandler.SendStr( FR_SURFACE   );
		break;
	  default:
		fSceneHandler.SendStr( FR_WIREFRAME );
		break;
	}

} 


//----- 
void G4FukuiRendererViewer::SendViewParameters () 
{
  // Calculates view representation based on extent of object being
  // viewed and (initial) direction of camera.  (Note: it can change
  // later due to user interaction via visualization system's GUI.)

#if defined DEBUG_FR_VIEW
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
    G4cout << "***** G4FukuiRendererViewer::SendViewParameters()" << G4endl;
#endif 

		//----- Magic number to decide camera distance automatically
	const    G4double        HOW_FAR            = 1000.0       ; // to define "infinity"
	const    G4double        MIN_HALF_ANGLE     = 0.01         ;
	const    G4double        MAX_HALF_ANGLE     = 0.499 * pi ;

		//----- (2A) CALC camera distance
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
		G4cout << "WARNING from FukuiRenderer (DAWN) driver:" << G4endl;
		G4cout << "  Camera cannot enter inside objects"      << G4endl;
	  }
		camera_distance = radius ; 
	}

		//----- (3A) CALC camera direction
	const G4Vector3D& camera_direction \
	  = fVP.GetViewpointDirection().unit();
	const G4double v_angle =  (180.0 / pi) * camera_direction.theta() ;
	const G4double h_angle =  (180.0 / pi) * camera_direction.phi  () ;

		//----- (2B), (3B) SEND camera position
	fSceneHandler.SendStrDouble3( FR_CAMERA_POSITION, 
			       camera_distance, 
			       v_angle, 
			       h_angle ); 

		//----- (4A) CALC target point
	const G4Point3D&  target_point
	  = fSceneHandler.GetScene()->GetStandardTargetPoint()
	  + fVP.GetCurrentTargetPoint();

		//----- (4B) SEND target point
	fSceneHandler.SendStrDouble3( FR_TARGET_POINT, 
			       target_point.x(), 
			       target_point.y(), 
			       target_point.z() );

		//----- (5A) CALC zoom factor
	const G4double   zoom_factor  = fVP.GetZoomFactor();

		//----- (5B) SEND zoom factor or focal length
	if( half_view_angle < MIN_HALF_ANGLE ) {

		const G4Point3D&  std_target_point \
	  		= fSceneHandler.GetScene()->GetStandardTargetPoint();

		fSceneHandler.SendStrDouble4( FR_ZOOM_FACTOR, 
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
		  = FR_HALF_SCREEN_SIZE / std::tan( half_view_angle ); 
		focal_distance *= zoom_factor ;
		fSceneHandler.SendStrDouble ( FR_FOCAL_DISTANCE, focal_distance );
	}

		//----- INVOKE GUI: not executed in the default setting
	if( fSceneHandler.GetSystem().IsGUIMode() ) {
			//----- send GUI command
		fSceneHandler.SendStr( FR_GUI );

			//----- wait the same command is sent back:
			//..... This avoids to send many data before
			//..... GUI session is over.
		fSceneHandler.GetPrimDest().WaitSendBack( FR_GUI );
	}

} 

#endif // G4VIS_BUILD_DAWN_DRIVER


