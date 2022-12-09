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
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#ifndef G4VTKVIEWER_HH
#define G4VTKVIEWER_HH

// #define G4VTKDEBUG  // Comment this out to suppress debug code.

#include "G4VViewer.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#include "vtkObject.h"
#include "vtkAutoInit.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleTerrain.h"
#include "vtkTextActor.h"
#pragma GCC diagnostic pop

VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType)

class vtkGeant4Callback : public vtkCommand
{
public:
  static vtkGeant4Callback* New() {return new vtkGeant4Callback;}

  vtkGeant4Callback() { fVP = nullptr; }
  void SetGeant4ViewParameters(G4ViewParameters *VP) { fVP = VP;}
  void SetVtkInitialValues(G4double parallelScaleIn, G4double cameraDistanceIn)
  {
    parallelScale = parallelScaleIn;
    cameraDistance = cameraDistanceIn;
  }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
  {
    vtkRenderer *ren = static_cast<vtkRenderer *>(caller);
    vtkCamera *cam = ren->GetActiveCamera();
    //G4cout << cam->GetFocalPoint()[0] << " " << cam->GetFocalPoint()[1] << " " << cam->GetFocalPoint()[2] << G4endl;
    //
    // G4cout << cam->GetPosition()[0] << " " << cam->GetPosition()[1] << " " << cam->GetPosition()[2] << G4endl;

    auto cp = cam->GetPosition();
    auto fp = cam->GetFocalPoint();
    auto ud = cam->GetViewUp();

    fVP->SetCurrentTargetPoint(G4Point3D(fp[0],fp[1],fp[2]));
    fVP->SetViewpointDirection((G4Point3D(cp[0],cp[1],cp[2]) -
                                G4Point3D(fp[0],fp[1],fp[2])).unit());
    fVP->SetUpVector(G4Vector3D(ud[0],ud[1], ud[2]));

    if(cam->GetParallelProjection()) {
      fVP->SetZoomFactor(parallelScale/cam->GetParallelScale());
    }
    else {
      auto cd = std::sqrt(std::pow(cp[0]-fp[0],2) +
                     std::pow(cp[1]-cp[1],2) +
                     std::pow(cp[2]-cp[2],2));
      fVP->SetZoomFactor(cameraDistance/cd);
    }
  }

protected:
  G4ViewParameters *fVP;
  G4double parallelScale;
  G4double cameraDistance;
};

class vtkInfoCallback : public vtkCommand
{
public:
  static vtkInfoCallback *New() { return new vtkInfoCallback; }

  vtkInfoCallback() {
    t1 = std::chrono::steady_clock::now();
    t2 = std::chrono::steady_clock::now();
  }
  void SetTextActor(vtkTextActor *txt) { this->TextActor = txt; }

  virtual void Execute(vtkObject *caller, unsigned long, void*)
  {
    vtkRenderer *ren = static_cast<vtkRenderer *>(caller);
    int      nActors = ren->GetActors()->GetNumberOfItems();
    vtkCamera   *cam = ren->GetActiveCamera();
    if(!cam) return;

    double      *pos     = cam->GetPosition();
    double      *foc     = cam->GetFocalPoint();
    double viewAngle     = cam->GetViewAngle();
    double distance      = cam->GetDistance();
    double parallelScale = cam->GetParallelScale();

    if(!pos) return;

    // Get current time
    t2 = std::chrono::steady_clock::now();

    // Frame rate calculation
    std::chrono::duration<double> tdiff = t2-t1;
    t1 = t2;
    float fps = 1.0/tdiff.count();

    // String for display
    snprintf(this->TextBuff,sizeof this->TextBuff,
	                   "camera position    : %.1f %.1f %.1f \n"
                           "camera focal point : %.1f %.1f %.1f \n"
                           "view angle         : %.1f\n"
                           "distance           : %.1f\n"
                           "parallel scale     : %.1f\n"
                           "number actors      : %i\n"
                           "fps                : %.1f",pos[0], pos[1], pos[2], foc[0], foc[1], foc[2], viewAngle, distance, parallelScale, nActors, fps);
    if(this->TextActor) {
      this->TextActor->SetInput(this->TextBuff);
    }
  }
protected:
  vtkTextActor *TextActor;
  char TextBuff[256];
  std::chrono::time_point<std::chrono::steady_clock> t1;
  std::chrono::time_point<std::chrono::steady_clock> t2;
};

class G4VtkViewer: public G4VViewer {
 public:
  G4VtkViewer(G4VSceneHandler&, const G4String& name);
  void Initialise();
  virtual ~G4VtkViewer();

  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();
  void FinishView();

  void DrawViewHUD();
  void DrawShadows();

  void ExportScreenShot(G4String, G4String);
  void ExportOBJScene(G4String);
  void ExportVRMLScene(G4String );
  void ExportVTPScene(G4String );

  void ExportView(){};
  void SetGeant4View(){};

  vtkNew<vtkTextActor> infoTextActor;
  vtkNew<vtkInfoCallback> infoCallback;
  vtkNew<vtkGeant4Callback> geant4Callback;
  vtkSmartPointer<vtkLight> light;
  vtkNew<vtkCamera> camera;
  vtkNew<vtkRenderer> renderer;
  vtkRenderWindow* _renderWindow;
  vtkRenderWindowInteractor* renderWindowInteractor;

 private:
  G4bool firstSetView = true;
};

#endif
