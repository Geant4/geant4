// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessCamera.cc,v 1.4 1999-12-15 14:54:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.
// Camera sub-menu.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"

void G4VisManMessenger::AddCommandCamera () {

  G4UIcommand* command;
  G4UIparameter* param;

  ///////////////////////////////////////  /vis~/camera/dolly  ////
  //camera \hline
  //camera /vis~/camera/dolly & distance &
  //camera Moves the camera in by this distance.
  //camera Reset with {\tt /vis~/camera/reset}. \\%
  command = new G4UIcommand ("/vis~/camera/dolly", this);
  command -> SetGuidance
    (
     "  Moves the camera in by this distance."
     "  Reset with /vis~/camera/reset."
     );
  param   =  new G4UIparameter ("move in", 'd', true);
  param   -> SetDefaultValue (0.);
  param   -> SetGuidance ("world coordinates");
  command -> SetParameter (param);
  fCommandList.append (command);

  ////////////////////////////////////////  /vis~/camera/orbit  ////
  //camera \hline
  //camera /vis~/camera/orbit & $N_\mathrm{frames}$, $\Delta\beta$ &
  //camera Orbits the scene about the up-vector, lights fixed to the scene.
  //camera Draws $N_\mathrm{frames}$ frames, the camera
  //camera rotated $\Delta\beta$ about the up-vector each frame. \\%
  command =  new G4UIcommand ("/vis~/camera/orbit", this);
  command -> SetGuidance
    (
     "Orbits the scene about the up-vector, lights fixed to the scene."
     "  Draws N frames, the camera "
     "rotated Delta-beta about the up-vector each frame."
     );
  param   =  new G4UIparameter ("No. of  frames", 'i', true);
  param   -> SetDefaultValue (1);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("angle_increment", 'd', true);
  param   -> SetGuidance ("degrees");
  param   -> SetDefaultValue (1.);
  command -> SetParameter (param);
  fCommandList.append (command);

  ///////////////////////////////////////  /vis~/camera/pan  ////
  //camera \hline
  //camera /vis~/camera/pan & right, up &
  //camera Moves the camera right and up by these amounts.
  //camera Reset with {\tt /vis~/camera/reset}. \\%
  command = new G4UIcommand ("/vis~/camera/pan", this);
  command -> SetGuidance
    (
     "  Moves the camera by right and up by these amounts."
     "  Reset with /vis~/camera/reset."
     );
  param   =  new G4UIparameter ("right shift", 'd', true);
  param   -> SetDefaultValue (0.);
  param   -> SetGuidance ("world coordinates");
  command -> SetParameter (param);
  param   =  new G4UIparameter ("up shift", 'd', true);
  param   -> SetDefaultValue (0.);
  param   -> SetGuidance ("world coordinates");
  command -> SetParameter (param);
  fCommandList.append (command);

  /////////////////////////////////////  /vis~/camera/projection_style  ////
  //camera \hline
  //camera /vis~/camera/projection\_style & choice &
  //camera Projection style (orthogonal, perspective). \\%
  command = new G4UIcommand ("/vis~/camera/projection_style", this);
  command -> SetGuidance
    (
     "Projection style (orthogonal, perspective)."
     );
  param   =  new G4UIparameter ("Style selector", 'i', true);
  param   -> SetDefaultValue  (-1);
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("Perspective view half angle.", 'd', true);
  param   -> SetDefaultValue  (0.);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ////////////////////////////////////////  /vis~/camera/spin  ////
  //camera \hline
  //camera /vis~/camera/spin & $N_\mathrm{frames}$, $\Delta\beta$ &
  //camera Spins the scene about the up-vector, lights fixed to the camera.
  //camera Draws $N_\mathrm{frames}$ frames, the scene
  //camera rotated $\Delta\beta$ about the up-vector each frame. \\%
  command =  new G4UIcommand ("/vis~/camera/spin", this);
  command -> SetGuidance
    (
     "Spins the scene about the up-vector, lights fixed to the camera."
     "  Draws N frames, the scene "
     "rotated Delta-beta about the up-vector each frame."
     );
  param   =  new G4UIparameter ("No. of  frames", 'i', true);
  param   -> SetDefaultValue (1);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("angle_increment", 'd', true);
  param   -> SetGuidance ("degrees");
  param   -> SetDefaultValue (1.);
  command -> SetParameter (param);
  fCommandList.append (command);

  ////////////////////////////////////////  /vis~/camera/viewpoint  ////
  //camera \hline
  //camera /vis~/camera/viewpoint & $\theta$ $\phi$ &
  //camera Set direction from target to camera as 
  //camera $\theta$, $\phi$ (in degrees).  Also changes lightpoint
  //camera direction if lights are set to move with camera.\\%
  command =  new G4UIcommand ("/vis~/camera/viewpoint", this);
  command -> SetGuidance
    (
     "Set direction from target to camera as "
     "theta, phi (in degrees).  Also changes lightpoint"
     "\n direction if lights are set to move with camera."
     );
  param   =  new G4UIparameter ("theta", 'd', true);
  param   -> SetGuidance ("degrees");
  param   -> SetDefaultValue (0.0);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("phi", 'd', true);
  param   -> SetGuidance ("degrees");
  param   -> SetDefaultValue (0.0);
  command -> SetParameter (param);
  fCommandList.append (command);

  ///////////////////////////////////////  /vis~/camera/window_size_hint  ////
  //camera \hline
  //camera /vis~/camera/window\_size\_hint & size&
  //camera Size in pixels (Square window `size $\Times$ size').  Affects
  //camera future create\_view commands (hint only). \\%
  command = new G4UIcommand    ("/vis~/camera/window_size_hint", this);
  command -> SetGuidance
    (
     "Size in pixels (Square window `size x size').  Affects "
     "future create_view commands (hint only)."
     );
  param   =  new G4UIparameter ("Window size hint in pixels", 'i', true);
  param   -> SetDefaultValue  (600);
  command -> SetParameter     (param);
  fCommandList.append (command);

  ///////////////////////////////////////  /vis~/camera/zoom  ////
  //camera \hline
  //camera /vis~/camera/zoom & factor &
  //camera Magnifies by this factor.
  //camera Reset with {\tt /vis~/camera/reset}. \\%
  command = new G4UIcommand ("/vis~/camera/zoom", this);
  command -> SetGuidance
    (
     "Magnifies by this factor."
     "  Reset with /vis~/camera/reset."
     );
  param   =  new G4UIparameter ("zoom factor", 'd', true);
  param   -> SetDefaultValue (1.);
  param   -> SetGuidance ("magnifies by this");
  command -> SetParameter (param);
  fCommandList.append (command);

}

void G4VisManMessenger::DoCommandCamera (const G4String& commandPath,
					 G4String& newValues) {

  ///////////////////////////////  /vis~/camera/dolly ////
  if (commandPath == "/vis~/camera/dolly") {
    if (ViewValid ()) {
      G4double in;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ; is >> in;
      G4cout << "Dolly " << in << " in." << G4endl;
      fpVMan -> SetCurrentViewParameters ().IncrementDolly (in);
      if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "Dolly distance changed to "
	     << fpVMan -> GetCurrentViewParameters ().GetDolly () << G4endl;
	if (fpVMan -> GetVerboseLevel () > 1) {
	  fpVMan -> PrintCurrentView ();
	}
      }
      // fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
      // fpVMan -> Draw ();
      // fpVMan -> Show ();
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (fpVMan -> GetCurrentViewParameters ());
	// Recalculate projection matrices, etc.
	pView -> SetView ();
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  ////////////////////////////////  /vis~/camera/orbit  ////
  if (commandPath  == "/vis~/camera/orbit") {
    if (ViewValid ()) {
      G4int nFrames;
      G4double dbeta;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ; is >> nFrames >> dbeta;
      dbeta = dbeta * deg;
      fpVMan -> SetCurrentViewParameters ().SetLightsMoveWithCamera (false);
      for (int i = 0; i < nFrames; i++) {
	RotateViewpointAboutUpVectorBy (dbeta);
	fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
	fpVMan -> Draw ();    
	fpVMan -> Show ();    
      }
    }
  }

  ///////////////////////////////  /vis~/camera/pan  ////
  if (commandPath == "/vis~/camera/pan") {
    if (ViewValid ()) {
      G4double right, up;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString);is >> right >> up;
      G4cout << "Pan " << right << " right, " << up << " up." << G4endl;
      fpVMan -> SetCurrentViewParameters ().Pan (right, up);
      if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "Current target point changed to "
	     << fpVMan -> GetCurrentViewParameters ().GetCurrentTargetPoint ()
	     << G4endl;
	if (fpVMan -> GetVerboseLevel () > 1) {
	  fpVMan -> PrintCurrentView ();
	}
      }
      // fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
      // fpVMan -> Draw ();    
      // fpVMan -> Show ();    
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (fpVMan -> GetCurrentViewParameters ());
	// Recalculate projection matrices, etc.
	pView -> SetView ();
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  /////////////////////////////////////  /vis~/camera/projection_style  ////
  if (commandPath == "/vis~/camera/projection_style") {
    G4int iStyle;
    G4double fieldHalfAngleDegrees;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> iStyle >> fieldHalfAngleDegrees;
    if (iStyle < 0 || iStyle > 1) {
      G4cout << "Available representation styles:";
      G4cout << "\n  0) orthogonal";
      G4cout << "\n  1) perspective (also specify field half angle in degrees)";
      G4cout << "\nChoose by specifying integer parameter and field half angle"
	" (if relevant).";
      G4cout << G4endl;
    }
    else {
      switch (iStyle) {
      default:
      case 0:
	fpVMan -> SetCurrentViewParameters ().SetFieldHalfAngle (0.);
	if (fpVMan -> GetVerboseLevel () > 0) {
	  G4cout << "\nOrthogonal projection set." << G4endl;
	}
	break;
      case 1:
	if (fieldHalfAngleDegrees == 0.0) {
	  fieldHalfAngleDegrees = 30.;
	  G4cout << "Field half angle defaulted to "
	       << fieldHalfAngleDegrees << " degrees." << G4endl;
	}
	if (fieldHalfAngleDegrees > 89.5 || fieldHalfAngleDegrees <= 0.0) {
	  G4cout << "Field half angle should be 0 < angle <= 89.5 degrees.";
	  G4cout << G4endl;
	}
	else {
	  fpVMan -> SetCurrentViewParameters ().
	    SetFieldHalfAngle (fieldHalfAngleDegrees * deg);
	  if (fpVMan -> GetVerboseLevel () > 0) {
	    G4cout << "\nPerspective projection set.";
	    G4cout << "\nField half angle set to " << fieldHalfAngleDegrees  
		 << " degrees" << " ("
		 << fpVMan -> GetCurrentViewParameters ().GetFieldHalfAngle ()
		 << " radians)." << G4endl;
	    if (fpVMan -> GetVerboseLevel () > 1) {
	      fpVMan -> PrintCurrentView ();
	    }
	  }
	  break;
	}
      }
      // if (ViewValid ()) {
	// fpVMan -> GetCurrentViewer () -> ClearView ();
	// fpVMan -> Draw ();    
	// fpVMan -> Show ();    
      // }
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (fpVMan -> GetCurrentViewParameters ());
	// Recalculate projection matrices, etc.
	pView -> SetView ();
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }

  ////////////////////////////////  /vis~/camera/spin  ////
  if (commandPath  == "/vis~/camera/spin") {
    if (ViewValid ()) {
      G4int nFrames;
      G4double dbeta;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ; is >> nFrames >> dbeta;
      dbeta = dbeta * deg;
      fpVMan -> SetCurrentViewParameters ().SetLightsMoveWithCamera (true);
      for (int i = 0; i < nFrames; i++) {
	RotateViewpointAboutUpVectorBy (dbeta);
	fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
	fpVMan -> Draw ();    
	fpVMan -> Show ();    
      }
    }
  }

  ////////////////////////////////////////  /vis~/camera/viewpoint  ////
  if (commandPath == "/vis~/camera/viewpoint") {
    G4double theta, phi ;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> theta >> phi;
    theta = theta * deg;
    phi   = phi   * deg;
    G4double x = sin (theta) * cos (phi);
    G4double y = sin (theta) * sin (phi);
    G4double z = cos (theta);
    G4Vector3D vp (x, y, z);
    fpVMan -> SetCurrentViewParameters().SetViewAndLights (vp);
    const G4ViewParameters& viewParams = fpVMan -> GetCurrentViewParameters ();
    if (fpVMan -> GetVerboseLevel () > 0) {
      G4cout << "Viewpoint direction set to " << vp << G4endl;
      if (viewParams.GetLightsMoveWithCamera ()) {
	G4cout << "Lightpoint direction set to "
	     << viewParams.GetLightpointDirection () << G4endl;
      }
      if (fpVMan -> GetVerboseLevel () > 1) {
	fpVMan -> PrintCurrentView ();
      }
    }
    // if (ViewValid ()) {
      // fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
      // fpVMan -> Draw ();    
      // fpVMan -> Show ();
    // }    
    G4VViewer* pView = fpVMan -> GetCurrentViewer ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (fpVMan -> GetCurrentViewParameters ());
      // Recalculate projection matrices, etc.
      pView -> SetView ();
    }
    G4cout << "Issue Draw or refresh to see effect." << G4endl;
  }

  ///////////////////////////////////////  /vis~/camera/window_size_hint  ////
  if (commandPath == "/vis~/camera/window_size_hint") {
    G4int size;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> size;
    fpVMan -> SetCurrentViewParameters ().SetWindowSizeHint (size, size);
    if (fpVMan -> GetVerboseLevel () > 0) {
      G4cout << "Window size hint set to " << size << " x " << size << G4endl;
      if (fpVMan -> GetVerboseLevel () > 1) {
	fpVMan -> PrintCurrentView ();
      }
    }
  }

  ////////////////////////////////  /vis~/camera/zoom  ////
  if (commandPath  == "/vis~/camera/zoom") {
    if (ViewValid ()) {
      G4double zoomBy;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString) ; is >> zoomBy;
      fpVMan -> SetCurrentViewParameters ().MultiplyZoomFactor (zoomBy);
      if (fpVMan -> GetVerboseLevel () > 0) {
	G4cout << "Zoom factor changed to "
	     << fpVMan -> GetCurrentViewParameters ().GetZoomFactor () << G4endl;
	if (fpVMan -> GetVerboseLevel () > 1) {
	  fpVMan -> PrintCurrentView ();
	}
      }
      // fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
      // fpVMan -> Draw ();    
      // fpVMan -> Show ();    
      G4VViewer* pView = fpVMan -> GetCurrentViewer ();
      if (pView) {
	// Copy current view parameters into current view.
	pView -> SetViewParameters (fpVMan -> GetCurrentViewParameters ());
	// Recalculate projection matrices, etc.
	pView -> SetView ();
      }
      G4cout << "Issue Draw or refresh to see effect." << G4endl;
    }
  }
}

void G4VisManMessenger::RotateViewpointAboutUpVectorBy (G4double dbeta) {
  // Rotates by fixed azimuthal angle dbeta.
  const G4ViewParameters& viewParams = fpVMan -> GetCurrentViewParameters ();
  G4Vector3D vp = viewParams.GetViewpointDirection ().unit ();
  G4Vector3D up = viewParams.GetUpVector ().unit ();
  G4Vector3D& zprime = up;
  G4double cosalpha = up.dot (vp);
  G4double sinalpha = sqrt (1. - pow (cosalpha, 2));
  G4Vector3D viewPoint = viewParams.GetViewpointDirection ().unit ();
  G4Vector3D yprime = (zprime.cross (viewPoint)).unit ();
  G4Vector3D xprime = yprime.cross (zprime);
  // Projection of vp on plane perpendicular to up...
  G4Vector3D a1 = sinalpha * xprime;
  // Required new projection...
  G4Vector3D a2 =
    sinalpha * (cos (dbeta) * xprime + sin (dbeta) * yprime);
  // Required Increment vector...
  G4Vector3D delta = a2 - a1;
  // So new viewpoint is...
  viewPoint += delta;
  fpVMan -> SetCurrentViewParameters ().SetViewAndLights (viewPoint);
}
