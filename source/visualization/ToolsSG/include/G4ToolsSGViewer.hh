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
// John Allison  6th October 2019

#if defined (G4VIS_BUILD_TOOLSSG_DRIVER) || defined (G4VIS_USE_TOOLSSG)

#ifndef G4TOOLSSGVIEWER_HH
#define G4TOOLSSGVIEWER_HH

#include "G4VViewer.hh"

#include "G4ToolsSGSceneHandler.hh"
#include "G4Scene.hh"

#include <tools/sg/device_interactor>
#include <tools/sg/separator>
#include <tools/sg/head_light>
#include <tools/sg/ortho>
#include <tools/sg/perspective>
#include <tools/sg/blend>
#include <tools/sg/noderef>
#include <tools/sg/keys>

template <class SG_SESSION,class SG_VIEWER>
class G4ToolsSGViewer : public G4VViewer, tools::sg::device_interactor {
  typedef G4VViewer parent;
  typedef tools::sg::device_interactor parent_interactor;
public: //tools::sg::device_interactor interface.
  virtual void key_press(const tools::sg::key_down_event& a_event) {
    fKeyPressed = true;
    fKeyShift = a_event.key() == tools::sg::key_shift()?true:false;
  }
  virtual void key_release(const tools::sg::key_up_event&) {fKeyPressed = false;}
  virtual void mouse_press(const tools::sg::mouse_down_event& a_event) {
    fMousePressed = true;
    fMousePressedX = a_event.x();
    fMousePressedY = a_event.y();
  }
  virtual void mouse_release(const tools::sg::mouse_up_event&) {fMousePressed = false;}    
  virtual void mouse_move(const tools::sg::mouse_move_event& a_event) {
    G4double x = a_event.x();
    G4double y = a_event.y();
    G4double dx = x-fMousePressedX;
    G4double dy = y-fMousePressedY;
    fMousePressedX = x;
    fMousePressedY = y;

    if (fMousePressed) {

      if (fKeyPressed && fKeyShift) {  // Translation (pan)

        const G4double sceneRadius = fSGSceneHandler.GetScene()->GetExtent().GetExtentRadius();
        const G4double scale = 300;  // Roughly pixels per window, empirically chosen
        const G4double dxScene = dx*sceneRadius/scale;
        const G4double dyScene = dy*sceneRadius/scale;
        fVP.IncrementPan(-dxScene,dyScene);

      } else {  // Rotation

        // Simple ad-hoc algorithms
        const G4Vector3D& x_prime = fVP.GetViewpointDirection().cross(fVP.GetUpVector());
        const G4Vector3D& y_prime = x_prime.cross(fVP.GetViewpointDirection());
        const G4double scale = 200;  // Roughly pixels per window, empirically chosen
        G4Vector3D newViewpointDirection = fVP.GetViewpointDirection();
        newViewpointDirection += dx*x_prime/scale;
        newViewpointDirection += dy*y_prime/scale;
        fVP.SetViewpointDirection(newViewpointDirection.unit());

        if (fVP.GetRotationStyle() == G4ViewParameters::freeRotation) {
          G4Vector3D newUpVector = fVP.GetUpVector();
          newUpVector += dx*x_prime/scale;
          newUpVector += dy*y_prime/scale;
          fVP.SetUpVector(newUpVector.unit());
        }
      }
    }

    SetView();
    DrawView();
  }
  virtual void wheel_rotate(const tools::sg::wheel_rotate_event& a_event) {
    const G4double angleY = a_event.angle();
    if (fVP.GetFieldHalfAngle() == 0.) {  // Orthographic projection
      const G4double scale = 500;  // Empirically chosen
      fVP.MultiplyZoomFactor(1.+angleY/scale);
    } else {                              // Perspective projection
      const G4double scale = fVP.GetFieldHalfAngle()/(10.*CLHEP::deg);  // Empirical
      fVP.SetDolly(fVP.GetDolly()+angleY/scale);
    }
    SetView();
    DrawView();
  }
public:
  G4ToolsSGViewer(SG_SESSION& a_session,G4ToolsSGSceneHandler& a_scene_handler, const G4String& a_name)
  :parent(a_scene_handler,a_scene_handler.IncrementViewCount(),a_name)
  ,fSGSession(a_session)
  ,fSGSceneHandler(a_scene_handler)
  ,fSGViewer(nullptr)
  ,fKeyPressed(false)
  ,fKeyShift(false)
  ,fMousePressed(false)
  ,fMousePressedX(0)
  ,fMousePressedY(0)
  {
    //::printf("debug : G4ToolsSGViewer::G4ToolsSGViewer.\n");
  }

  virtual ~G4ToolsSGViewer() {
    delete fSGViewer;
    //::printf("debug : G4ToolsSGViewer::~G4ToolsSGViewer.\n");
  }
protected:
  G4ToolsSGViewer(const G4ToolsSGViewer& a_from)
  :parent(a_from)
  ,parent_interactor(a_from)
  ,fSGSession(a_from.fSGSession)
  ,fSGSceneHandler(a_from.fSGSceneHandler)
  ,fSGViewer(nullptr)
  ,fKeyPressed(false)
  ,fKeyShift(false)
  ,fMousePressed(false)
  ,fMousePressedX(0)
  ,fMousePressedY(0)
  {}
  G4ToolsSGViewer& operator=(const G4ToolsSGViewer&) {return *this;}
public:  
  virtual void Initialise() {
    if(fSGViewer) return; //done.
    fVP.SetAutoRefresh(true);
    fDefaultVP.SetAutoRefresh(true);
    //::printf("debug : G4ToolsSGViewer::Initialise\n");
    //////////////////////////////////////////////////////////
    /// create the viewer, set the scene graph ///////////////
    //////////////////////////////////////////////////////////
    fSGViewer = new SG_VIEWER(fSGSession
      ,fVP.GetWindowAbsoluteLocationHintX(1440)
      ,fVP.GetWindowAbsoluteLocationHintY(900)
      ,fVP.GetWindowSizeHintX()
      ,fVP.GetWindowSizeHintY()
      ,fName);
    if(!fSGViewer->has_window()) {
      fViewId = -1;  // This flags an error.
      G4cerr << "G4ToolsSGViewer::Initialise : SG_VIEWER::has_window() failed." << G4endl;
      return;
    }

    fSGViewer->set_clear_color(0.3,0.3,0.3,1);
    fSGViewer->set_device_interactor(this);
  }
  
  virtual void SetView() {
    //::printf("debug : G4ToolsSGViewer::SetView\n");
    if(!fSceneHandler.GetScene()) {
      fSGViewer->set_clear_color(0.3,0.3,0.3,1); //some grey color to signal the user that something is wrong.
      G4cerr << "G4ToolsSGViewer::SetView : no G4Scene.." << G4endl;
      return;
    }
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // Get radius of scene, etc.
    // Note that this procedure properly takes into account zoom, dolly and pan.
    const G4Point3D targetPoint
      = fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint ();
    G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
    if(radius<=0.) radius = 1.;
    const G4double cameraDistance = fVP.GetCameraDistance (radius);
    const G4Point3D cameraPosition = targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
    const G4Normal3D& up = fVP.GetUpVector ();  
    const G4double pnear  = fVP.GetNearDistance (cameraDistance, radius);
    const G4double pfar   = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  //const G4double right  = fVP.GetFrontHalfHeight (pnear, radius);
  //const G4double left   = -right;
    const G4double top    = fVP.GetFrontHalfHeight (pnear, radius);
    const G4double bottom = -top;
    // sanity check :
    tools::vec3f dir(float(targetPoint.x()-cameraPosition.x()),
                     float(targetPoint.y()-cameraPosition.y()),
                     float(targetPoint.z()-cameraPosition.z()));
    if(!dir.length()) {
      fSGViewer->set_clear_color(0.3,0.3,0.3,1);
      G4cerr << "G4ToolsSGViewer::SetView : null size viewer area." << G4endl;
      return;      
    }
    //::printf("debug : GB : SetView : 1-001 : %g %g %g : %g %g %g : h %g : %g %g\n",
    //  targetPoint.x(),targetPoint.y(),targetPoint.z(),
    //  cameraPosition.x(),cameraPosition.y(),cameraPosition.z(),top-bottom,pnear,pfar);

    //////////////////////////////////////////////////////////
    /// create scene graph ///////////////////////////////////
    //////////////////////////////////////////////////////////
    // Set projection, then create the tools::sg camera node :
    tools::sg::base_camera* _camera = nullptr;
    if (fVP.GetFieldHalfAngle() <= 0.) {
      if((top-bottom)<=0) {
        fSGViewer->set_clear_color(0.3,0.3,0.3,1);
        G4cerr << "G4ToolsSGViewer::SetView : for ortho camera, (top-bottom)<=0." << G4endl;
        return;
      }
      tools::sg::ortho* ortho_camera = new tools::sg::ortho;
      ortho_camera->height.value(float(top-bottom));
      _camera = ortho_camera;
    } else {
      tools::sg::perspective* perspective_camera = new tools::sg::perspective;
      perspective_camera->height_angle.value(float(2*fVP.GetFieldHalfAngle()));
      _camera = perspective_camera;
    }
    
    _camera->position.value
      (tools::vec3f(float(cameraPosition.x()),
		    float(cameraPosition.y()),
		    float(cameraPosition.z())));
    _camera->znear.value(float(pnear));
    _camera->zfar.value(float(pfar));

    _camera->look_at(dir,tools::vec3f(up.x(),up.y(),up.z()));  //same logic as in G4OpenInventorViewer.
    
    tools::sg::separator* sep = new tools::sg::separator;
    sep->add(_camera);

   {tools::sg::head_light* light = new tools::sg::head_light;
  //light->direction = tools::vec3f(1,-1,-10);
    light->direction = tools::vec3f(0,0,-1);
    light->on = true;
    sep->add(light);}

   {tools::sg::blend* blend = new tools::sg::blend;
    blend->on = true; //to handle transparency.
    sep->add(blend);}
 
    sep->add(new tools::sg::noderef(*fSGSceneHandler.ToolsSGScene()));
    
    fSGViewer->sg().clear();
    fSGViewer->sg().add(sep); //give sep ownership to the viewer.
    
   {G4Color background = fVP.GetBackgroundColour ();
    fSGViewer->set_clear_color(float(background.GetRed()),float(background.GetBlue()),float(background.GetBlue()),1);}
  }

  virtual void ClearView() {}

  virtual void DrawView() {
    if (!fNeedKernelVisit) KernelVisitDecision();
    G4bool kernelVisitWasNeeded = fNeedKernelVisit; // Keep (ProcessView resets).
    fLastVP = fVP;
    ProcessView();  // Clears store and processes scene only if necessary.
    if (kernelVisitWasNeeded) {
      // We might need to do something if the kernel was visited.
    } else {
    }
    FinishView ();       // Flush streams and/or swap buffers.
  }

  virtual void ShowView() {
    if(!fSGViewer) return;
    fSGViewer->show();
    fSGViewer->win_render();
    fSGSession.sync();
  }

  virtual void FinishView() {
    if(!fSGViewer) return;
    fSGViewer->show();
    fSGViewer->win_render();
    fSGSession.sync();
  }


#ifdef G4MULTITHREADED
  virtual void SwitchToVisSubThread() {}
  
  virtual void SwitchToMasterThread() {
    if (G4Threading::IsMultithreadedApplication()) {
      // I have not figured out how to draw during a run.
      //
      // Setting fNeedKernelVisit=true causes scene deletion and a complete rebuild,
      // including trajectories, hits, etc. from kept events.
      //
      // Clearly this is a limitation because even if you run 1000 events you only
      // get those kept (default 100), and even worse, if end-if-event-action is
      // "refresh", you only get one event (the last I think).
      //
      // Also, strictly, there is no need to rebuid run-duration models (detector),
      // but a complete rebuild is the easiest way (already imeplemented).
      fNeedKernelVisit = true;
      DrawView();  // Draw trajectories, etc., from kept events
    }
  }
#endif
  
  //SG_VIEWER* sg_viewer() {return fSGViewer;}
protected:
  void KernelVisitDecision () {
    if (CompareForKernelVisit(fLastVP)) {
      NeedKernelVisit ();  // Sets fNeedKernelVisit.
    }
  }
  
  G4bool CompareForKernelVisit(G4ViewParameters& lastVP) {
    // Typical comparison.  Taken from OpenGL.
    if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.GetNumberOfCloudPoints()  != fVP.GetNumberOfCloudPoints())  ||
      (lastVP.IsAuxEdgeVisible ()   != fVP.IsAuxEdgeVisible ())   ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (lastVP.GetCBDAlgorithmNumber() !=
       fVP.GetCBDAlgorithmNumber())                               ||
      (lastVP.IsSection ()          != fVP.IsSection ())          ||
      (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
      (lastVP.GetGlobalMarkerScale()    != fVP.GetGlobalMarkerScale())    ||
      (lastVP.GetGlobalLineWidthScale() != fVP.GetGlobalLineWidthScale()) ||
      (lastVP.IsMarkerNotHidden ()  != fVP.IsMarkerNotHidden ())  ||
      (lastVP.GetDefaultVisAttributes()->GetColour() !=
       fVP.GetDefaultVisAttributes()->GetColour())                ||
      (lastVP.GetDefaultTextVisAttributes()->GetColour() !=
       fVP.GetDefaultTextVisAttributes()->GetColour())            ||
      (lastVP.GetBackgroundColour ()!= fVP.GetBackgroundColour ())||
      (lastVP.IsPicking ()          != fVP.IsPicking ())          ||
      (lastVP.GetVisAttributesModifiers() !=
       fVP.GetVisAttributesModifiers())                           ||
      (lastVP.IsSpecialMeshRendering() !=
       fVP.IsSpecialMeshRendering())
      ) {
      return true;
    }

    if (lastVP.IsDensityCulling () &&
       (lastVP.GetVisibleDensity () != fVP.GetVisibleDensity ())) return true;

    if (lastVP.GetCBDAlgorithmNumber() > 0) {
      if (lastVP.GetCBDParameters().size() != fVP.GetCBDParameters().size()) return true;
      else if (lastVP.GetCBDParameters() != fVP.GetCBDParameters()) return true;
    }

    if (lastVP.IsExplode () &&
       (lastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
      return true;

    if (lastVP.IsSpecialMeshRendering() &&
	(lastVP.GetSpecialMeshVolumes() != fVP.GetSpecialMeshVolumes()))
      return true;

    return false;
  }

//  void keyPressEvent        (KeyEvent*);
//  void keyReleaseEvent      (KeyEvent*);
//  void mouseDoubleClickEvent(MouseEvent*);
//  void mouseMoveEvent       (MouseEvent*);
//  void mousePressEvent      (MouseEvent*);
//  void mouseReleaseEvent    (MouseEvent*);
//  void wheelEvent           (WheelEvent*);

protected:
  SG_SESSION& fSGSession;
  G4ToolsSGSceneHandler& fSGSceneHandler;
  SG_VIEWER* fSGViewer;
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
  
  G4bool fKeyPressed;
  G4bool fKeyShift;
  G4bool fMousePressed;
  G4double fMousePressedX, fMousePressedY;
};

#endif

#endif
