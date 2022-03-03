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

#ifndef G4TOOLSSGVIEWER_HH
#define G4TOOLSSGVIEWER_HH

#include "G4VViewer.hh"

#include "G4ToolsSGSceneHandler.hh"
#include "G4Scene.hh"
#include "G4VVisCommand.hh"

#include <tools/sg/device_interactor>
#include <tools/sg/separator>
#include <tools/sg/ortho>
#include <tools/sg/perspective>
#include <tools/sg/torche>
#include <tools/sg/blend>
#include <tools/sg/noderef>
#include <tools/sg/keys>

#include <tools/tokenize>
#include <tools/sg/write_paper>

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
    //::printf("debug : G4ToolsSGViewer::G4ToolsSGViewer: %lu, %s\n",this,a_name.c_str());
    Messenger::Create();
  }

  virtual ~G4ToolsSGViewer() {
    //::printf("debug : G4ToolsSGViewer::~G4ToolsSGViewer: %lu\n",this);
    //WARNING : nodes may refer f_gl2ps_mgr, f_zb_mgr (to handle gstos (for GPU) or textures), then
    //          we have to delete them first.
    fSGViewer->sg().clear();
    delete fSGViewer;
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
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    /*
    G4cout << "debug : 0002 : radius " << radius << std::endl;
    G4cout << "debug : cameraDistance : " << cameraDistance << std::endl;
    G4cout << "debug : fieldHalfAngle : " << fVP.GetFieldHalfAngle() << std::endl;
    G4cout << "debug : zoomFactor : " << fVP.GetZoomFactor() << std::endl;
    G4cout << "debug : up : " << up.x() << " " << up.y() << " " << up.z() << std::endl;
    G4cout << "debug : targetPoint : " << targetPoint.x() << " " << targetPoint.y() << " " << targetPoint.z() << std::endl;
    G4cout << "debug : cameraPosition : " << cameraPosition.x() << " " << cameraPosition.y() << " " << cameraPosition.z() << std::endl;
    G4cout << "debug : camera : znear " << pnear << ", zfar " << pfar << std::endl;
    */
    //////////////////////////////////////////////////////////
    /// create scene graph ///////////////////////////////////
    //////////////////////////////////////////////////////////
    // Set projection, then create the tools::sg camera node :
    tools::sg::base_camera* _camera = nullptr;
    if (fVP.GetFieldHalfAngle() <= 0.) {
      //G4cout << "debug : camera : ortho : top " << top << " bottom " << bottom << " top-bottom " << top-bottom << std::endl;
      if((top-bottom)<=0) {
        fSGViewer->set_clear_color(0.3,0.3,0.3,1);
        G4cerr << "G4ToolsSGViewer::SetView : for ortho camera, (top-bottom)<=0." << G4endl;
        return;
      }
      tools::sg::ortho* ortho_camera = new tools::sg::ortho;
      ortho_camera->height.value(float(top-bottom));
      _camera = ortho_camera;
    } else {
      //G4cout << "debug : camera : perspec : heightAngle " << float(2*fVP.GetFieldHalfAngle()) << std::endl;
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

    /*
    const G4Vector3D& lightDirection = fVP.GetLightpointDirection();
    G4cout << "debug : lightDirection : " << lightDirection.x() << " " << lightDirection.y() << " " << lightDirection.z() << std::endl;
    const G4Vector3D& actualLightDirection = fVP.GetActualLightpointDirection();
    G4cout << "debug : actualLightDirection : " << actualLightDirection.x() << " " << actualLightDirection.y() << " " << actualLightDirection.z() << std::endl;
    */

    CreateSG(_camera,fVP.GetActualLightpointDirection());
    
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

  virtual void ShowView() {FinishView();}

  virtual void FinishView() {
    if(fSGViewer) {
      fSGSceneHandler.TouchPlotters(fSGViewer->sg());
      fSGViewer->show();
      fSGViewer->win_render();
      fSGSession.sync();
    }
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
  void CreateSG(tools::sg::base_camera* a_camera,const G4Vector3D& a_light_dir) {
    tools::sg::group& _parent = fSGViewer->sg();
    _parent.clear();    

    ///////////////////////////////////////////////////
    /// 2D scene graph: ///////////////////////////////
    ///////////////////////////////////////////////////
    tools::sg::separator* scene_2D = new tools::sg::separator;
    _parent.add(scene_2D);
    scene_2D->add(new tools::sg::noderef(fSGSceneHandler.GetTransient2DObjects()));
    scene_2D->add(new tools::sg::noderef(fSGSceneHandler.GetPersistent2DObjects()));
  
    ///////////////////////////////////////////////////
    /// 3D scene graph: ///////////////////////////////
    ///////////////////////////////////////////////////
    tools::sg::separator* scene_3D = new tools::sg::separator;
    _parent.add(scene_3D);
  
    scene_3D->add(a_camera);
  
   {tools::sg::torche* light = new tools::sg::torche;
    light->on = true;
    light->direction = tools::vec3f(-a_light_dir.x(),-a_light_dir.y(),-a_light_dir.z());
    scene_3D->add(light);}
  
   {tools::sg::blend* blend = new tools::sg::blend;
    blend->on = true; //to handle transparency.
    scene_3D->add(blend);}

    scene_3D->add(new tools::sg::noderef(fSGSceneHandler.GetTransient3DObjects()));
    scene_3D->add(new tools::sg::noderef(fSGSceneHandler.GetPersistent3DObjects()));
  }
  
  void Export(const G4String& a_format,const G4String& a_file) {
    if(!fSGViewer) return;
    const G4Colour& back_color = fVP.GetBackgroundColour();
    if(!write_paper(G4cout,f_gl2ps_mgr,f_zb_mgr,0,0,
                    float(back_color.GetRed()),float(back_color.GetGreen()),float(back_color.GetBlue()),float(back_color.GetAlpha()),
		    fSGViewer->sg(),
		    fSGViewer->width(),fSGViewer->height(),a_file,a_format)) {
      G4cout << "G4ToolsSGViewer::Export: write_paper() failed." << G4endl;
      return;
    }
  }

protected:  
  class Messenger: public G4VVisCommand {
  public:  
    static void Create() {static Messenger s_messenger;}
  private:  
    Messenger() {
      G4UIparameter* parameter;
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      write_scene = new G4UIcommand("/vis/tsg/export", this);
      write_scene->SetGuidance("Write the content of the current viewer in a file at various formats.");
      write_scene->SetGuidance("Default file is out.eps and default format is gl2ps_eps.");
      write_scene->SetGuidance("Available formats are:");
      write_scene->SetGuidance("- gl2ps_eps: gl2ps producing eps");
      write_scene->SetGuidance("- gl2ps_ps:  gl2ps producing ps");
      write_scene->SetGuidance("- gl2ps_pdf: gl2ps producing pdf");
      write_scene->SetGuidance("- gl2ps_svg: gl2ps producing svg");
      write_scene->SetGuidance("- gl2ps_tex: gl2ps producing tex");
      write_scene->SetGuidance("- gl2ps_pgf: gl2ps producing pgf");
      write_scene->SetGuidance("- zb_ps: tools::sg offscreen zbuffer put in a PostScript file.");

      parameter = new G4UIparameter("format",'s',true);
      parameter->SetDefaultValue("gl2ps_eps");
      write_scene->SetParameter (parameter);

      parameter = new G4UIparameter("file",'s',true);
      parameter->SetDefaultValue("out.eps");
      write_scene->SetParameter (parameter);
    }
    virtual ~Messenger() {
      delete write_scene;
    }
  public:
    virtual void SetNewValue(G4UIcommand* a_cmd,G4String a_value) {
      G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
      G4VViewer* viewer = fpVisManager -> GetCurrentViewer ();
      if (!viewer) {
        if (verbosity >= G4VisManager::errors) G4cerr << "ERROR: No current viewer." << G4endl;
        return;
      }
      G4ToolsSGViewer* tsg_viewer = dynamic_cast<G4ToolsSGViewer*>(viewer);
      if(!tsg_viewer) {
        G4cout << "G4ToolsSGViewer::SetNewValue:"
               << " current viewer is not a G4ToolsSGViewer." << G4endl;
        return;
      }
      std::vector<std::string> args;
      tools::double_quotes_tokenize(a_value,args);
      if(args.size()!=a_cmd->GetParameterEntries()) return;
      if(a_cmd==write_scene) {
        tsg_viewer->Export(args[0],args[1]);
      }
    }
  private:
    G4UIcommand* write_scene;
  };
  
protected:
  SG_SESSION& fSGSession;
  G4ToolsSGSceneHandler& fSGSceneHandler;
  SG_VIEWER* fSGViewer;
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
  
  G4bool fKeyPressed;
  G4bool fKeyShift;
  G4bool fMousePressed;
  G4double fMousePressedX, fMousePressedY;

  tools::sg::zb_manager f_zb_mgr;
  tools::sg::gl2ps_manager f_gl2ps_mgr;
  
};

#endif
