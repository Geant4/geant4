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

#include "G4VtkViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4VtkSceneHandler.hh"

#include "vtkRendererCollection.h"
#include "vtkLightCollection.h"

#include "vtkWindowToImageFilter.h"
#include "vtkImageWriter.h"
#include "vtkBMPWriter.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"
#include "vtkPNMWriter.h"
#include "vtkTIFFWriter.h"
#include "vtkPostScriptWriter.h"
#include "vtkOBJExporter.h"
#include "vtkVRMLExporter.h"
#include "vtkSingleVTPExporter.h"

#include "vtkShadowMapPass.h"
#include "vtkShadowMapBakerPass.h"
#include "vtkSequencePass.h"
#include "vtkCameraPass.h"
#include "vtkRenderPass.h"
#include "vtkRenderPassCollection.h"

#include "vtkOpenGLRenderer.h"

G4VtkViewer::G4VtkViewer(G4VSceneHandler& sceneHandler, const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
{
  vtkObject::GlobalWarningDisplayOff();

  // Set default and current view parameters
  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);
}

void G4VtkViewer::Initialise()
{
  _renderWindow          = vtkRenderWindow::New();
  renderWindowInteractor = vtkRenderWindowInteractor::New();

#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::G4VtkViewer" << G4endl;
  G4cout << "G4VtkViewer::G4VtkViewer> " << fVP.GetWindowSizeHintX() << " "
  << fVP.GetWindowSizeHintY() << G4endl;
  G4cout << "G4VtkViewer::G4VtkViewer> " << fVP.GetWindowLocationHintX() << " "
  << fVP.GetWindowLocationHintY() << G4endl;
#endif

  // Need windowSizeX/Y - obtain from _renderWindow?
  G4int screenSizeX = _renderWindow->GetScreenSize()[0];
  G4int screenSizeY = _renderWindow->GetScreenSize()[1];
  G4int positionX = fVP.GetWindowLocationHintX();
  if (fVP.IsWindowLocationHintXNegative()) {
    positionX = screenSizeX + positionX - fVP.GetWindowSizeHintX();
  }
  G4int positionY = fVP.GetWindowLocationHintY();
  if (!fVP.IsWindowLocationHintYNegative()) {
    positionY = screenSizeY + positionY - fVP.GetWindowSizeHintY();
  }
  _renderWindow->SetPosition(positionX, positionY);
#ifdef __APPLE__
  // Adjust window size for Apple to make it correspond to OpenGL.
  // Maybe it's OpenGL that shoud be adjusted.
  const G4double pixelFactor = 2.;
#else
  const G4double pixelFactor = 1.;
#endif
  _renderWindow->SetSize
  (pixelFactor*fVP.GetWindowSizeHintX(),pixelFactor*fVP.GetWindowSizeHintY());
  _renderWindow->SetWindowName("Vtk viewer");

  _renderWindow->AddRenderer(renderer);
  renderWindowInteractor->SetRenderWindow(_renderWindow);

  // TODO proper camera parameter settings
  camera->SetPosition(0, 0, 1000);
  camera->SetFocalPoint(0, 0, 0);
  renderer->SetActiveCamera(camera);

  //renderer->SetUseHiddenLineRemoval(1);  // TODO needs to be an option
  //renderer->SetUseShadows(1);            // TODO needs to be an option

  // Set callback to match VTK parameters to Geant4
  geant4Callback->SetGeant4ViewParameters(&fVP);
  renderer->AddObserver(vtkCommand::EndEvent, geant4Callback);

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);

  // DrawShadows();
}

G4VtkViewer::~G4VtkViewer() {}

void G4VtkViewer::SetView() {

  // background colour
  const G4Colour backgroundColour = fVP.GetBackgroundColour();
  renderer->SetBackground(backgroundColour.GetRed(), backgroundColour.GetGreen(), backgroundColour.GetBlue());

  // target and camera positions
  G4double radius = fSceneHandler.GetExtent().GetExtentRadius();
  if(radius <= 0.)
    {radius = 1.;}
  G4double cameraDistance      = fVP.GetCameraDistance(radius);
  G4Point3D viewpointDirection = fVP.GetViewpointDirection();
  G4Point3D targetPoint        = fVP.GetCurrentTargetPoint();
  G4Point3D cameraPosition =
    targetPoint + viewpointDirection.unit() * cameraDistance;
  renderer->GetActiveCamera()->SetFocalPoint(targetPoint.x(),
                                             targetPoint.y(),
                                             targetPoint.z());
  renderer->GetActiveCamera()->SetPosition(cameraPosition.x(),
                                           cameraPosition.y(),
                                           cameraPosition.z());
  renderer->GetActiveCamera()->SetParallelScale(cameraDistance);

  // need to set camera distance and parallel scale on first set view
  if(firstSetView)
  {
    geant4Callback->SetVtkInitialValues(cameraDistance, cameraDistance);
    firstSetView = false;
  }

  // projection type and view angle and zoom factor
  G4double      fieldHalfAngle = fVP.GetFieldHalfAngle();
  G4double          zoomFactor = fVP.GetZoomFactor();
  vtkCamera* activeCamera = renderer->GetActiveCamera();
  if(fieldHalfAngle == 0) {
    activeCamera->SetParallelProjection(1);
    activeCamera->SetParallelScale(activeCamera->GetParallelScale()/zoomFactor);
  }
  else {
    activeCamera->SetParallelProjection(0);
    activeCamera->SetViewAngle(2*fieldHalfAngle/M_PI*180);
    activeCamera->SetPosition(cameraPosition.x()/zoomFactor,
                                             cameraPosition.y()/zoomFactor,
                                             cameraPosition.z()/zoomFactor);
  }

  // camera roll
  // renderer->GetActiveCamera()->SetRoll(0);

  // camera up direction
  const G4Vector3D upVector = fVP.GetUpVector();
  renderer->GetActiveCamera()->SetViewUp(upVector.x(),
                                         upVector.y(),
                                         upVector.z());

  // Light
  const G4Vector3D lightDirection = fVP.GetLightpointDirection();
  G4bool lightsMoveWithCamera     = fVP.GetLightsMoveWithCamera();
  G4Vector3D lightPosition =
    targetPoint + lightDirection.unit() * cameraDistance;

  vtkLightCollection* currentLights = renderer->GetLights();
  if (currentLights->GetNumberOfItems() != 0)
  {
    auto currentLight = dynamic_cast<vtkLight*>(currentLights->GetItemAsObject(0));
    if (currentLight)
    {
      currentLight->SetPosition(lightPosition.x(),
                                lightPosition.y(),
                                lightPosition.z());
      if (lightsMoveWithCamera)
      {currentLight->SetLightTypeToCameraLight();}
      else
      {currentLight->SetLightTypeToSceneLight();}
    }
  }

  // Rotation style
#if 0
  G4ViewParameters::RotationStyle rotationStyle  = fVP.GetRotationStyle();
  if (rotationStyle == G4ViewParameters::RotationStyle::freeRotation) {
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    _renderWindow->GetInteractor()->SetInteractorStyle(style);
  }
  else if(rotationStyle == G4ViewParameters::RotationStyle::constrainUpDirection) {
    // camera->SetViewUp(upVector.x(), upVector.y(), upVector.z());
    vtkSmartPointer<vtkInteractorStyleTerrain> style =
      vtkSmartPointer<vtkInteractorStyleTerrain>::New();
    _renderWindow->GetInteractor()->SetInteractorStyle(style);
  }
#endif
}

void G4VtkViewer::ClearView() {
  vtkActorCollection *actors = renderer->GetActors();
  vtkActor *actor = actors->GetLastActor();

  while(actor) {
#ifdef G4VTKDEBUG
    G4cout << "G4VtkViewer::ClearView() remove actor " << actor << G4endl;
#endif
    renderer->RemoveActor(actor);
    actor = actors->GetLastActor();
  }

  vtkPropCollection *props = renderer->GetViewProps();
  vtkProp *prop  = props->GetLastProp();

  while(prop) {
#ifdef G4VTKDEBUG
    G4cout << "G4VtkViewer::ClearView() remove prop " << prop << G4endl;
#endif
    renderer->RemoveViewProp(prop);
    prop = props->GetLastProp();
  }
}

void G4VtkViewer::DrawView() {
  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  NeedKernelVisit();  // Default is - always visit G4 kernel.
  // Note: this routine sets the fNeedKernelVisit flag of *all* the
  // views of the scene.

  ProcessView();  // The basic logic is here.

  // Add HUD
  DrawViewHUD();

  // ...before finally...
  FinishView();  // Flush streams and/or swap buffers.
}

void G4VtkViewer::DrawViewHUD()
{
  // make sure text is always visible
  G4Colour colour = fVP.GetBackgroundColour();
  infoTextActor->GetTextProperty()->SetColor(std::fmod(colour.GetRed() + 0.5, 1.0),
                                             std::fmod(colour.GetGreen() + 0.5, 1.0),
                                             std::fmod(colour.GetBlue() + 0.5, 1.0));
  infoTextActor->GetTextProperty()->SetFontSize(20);
  infoCallback->SetTextActor(infoTextActor);
  renderer->AddObserver(vtkCommand::EndEvent, infoCallback);
  renderer->AddActor(infoTextActor);
}

void G4VtkViewer::DrawShadows()
{
  _renderWindow->SetMultiSamples(0);

  vtkNew<vtkShadowMapPass> shadows;
  vtkNew<vtkSequencePass> seq;

  vtkNew<vtkRenderPassCollection> passes;
  passes->AddItem(shadows->GetShadowMapBakerPass());
  passes->AddItem(shadows);
  seq->SetPasses(passes);

  vtkNew<vtkCameraPass> cameraP;
  cameraP->SetDelegatePass(seq);

  // tell the renderer to use our render pass pipeline
  vtkOpenGLRenderer* glrenderer = dynamic_cast<vtkOpenGLRenderer*>(renderer.GetPointer());
  glrenderer->SetPass(cameraP);
}

void G4VtkViewer::ShowView()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ShowView() called." << G4endl;
  // static_cast<G4VtkSceneHandler&>(fSceneHandler).PrintStores();
#endif

  G4VtkSceneHandler& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Modified();

  infoTextActor->GetTextProperty()->SetFontSize(28);
  G4Colour colour = fVP.GetBackgroundColour();

  // make sure text is always visible
  infoTextActor->GetTextProperty()->SetColor(std::fmod(colour.GetRed() + 0.5, 1.0),
                                             std::fmod(colour.GetGreen() + 0.5, 1.0),
                                             std::fmod(colour.GetBlue() + 0.5, 1.0));
  infoTextActor->GetTextProperty()->SetFontSize(20);
  infoCallback->SetTextActor(infoTextActor);
  renderer->AddObserver(vtkCommand::EndEvent, infoCallback);
  geant4Callback->SetGeant4ViewParameters(&fVP);
  renderer->AddObserver(vtkCommand::EndEvent, geant4Callback);
  renderer->AddActor(infoTextActor);
}

void G4VtkViewer::FinishView()
{
  G4VtkSceneHandler& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Modified();

  _renderWindow->Render();
  _renderWindow->GetInteractor()->Initialize();
  _renderWindow->GetInteractor()->Start();
}

void G4VtkViewer::ExportScreenShot(G4String path, G4String format)
{

  vtkImageWriter *imWriter = nullptr;

  if(format == "bmp") {
    imWriter = vtkBMPWriter::New();
  }
  else if (format == "jpg") {
    imWriter = vtkJPEGWriter::New();
  }
  else if (format == "pnm") {
    imWriter = vtkPNMWriter::New();
  }
  else if (format == "png") {
    imWriter = vtkPNGWriter::New();
  }
  else if (format == "tiff") {
    imWriter = vtkTIFFWriter::New();
  }
  else if (format == "ps") {
    imWriter = vtkPostScriptWriter::New();
  }
  else {
    imWriter = vtkPNGWriter::New();
  }

  _renderWindow->Render();

  vtkSmartPointer<vtkWindowToImageFilter> winToImage = vtkSmartPointer<vtkWindowToImageFilter>::New();
  winToImage->SetInput(_renderWindow);
  winToImage->SetScale(1);
  if(format == "ps")
  {
    winToImage->SetInputBufferTypeToRGB();
    winToImage->ReadFrontBufferOff();
    winToImage->Update();
  }
  else
  {winToImage->SetInputBufferTypeToRGBA();}

  imWriter->SetFileName((path+"."+format).c_str());
  imWriter->SetInputConnection(winToImage->GetOutputPort());
  imWriter->Write();
}

void G4VtkViewer::ExportOBJScene(G4String path)
{
  vtkSmartPointer<vtkRenderWindow> _rw1 = vtkSmartPointer<vtkRenderWindow>::New();
  _rw1->AddRenderer(_renderWindow->GetRenderers()->GetFirstRenderer());
  vtkSmartPointer<vtkOBJExporter> exporter = vtkSmartPointer<vtkOBJExporter>::New();
  exporter->SetRenderWindow(_rw1);
  exporter->SetFilePrefix(path.c_str());
  exporter->Write();
}

void G4VtkViewer::ExportVRMLScene(G4String path)
{
  vtkSmartPointer<vtkRenderWindow> _rw1 = vtkSmartPointer<vtkRenderWindow>::New();
  _rw1->AddRenderer(_renderWindow->GetRenderers()->GetFirstRenderer());
  vtkSmartPointer<vtkVRMLExporter> exporter = vtkSmartPointer<vtkVRMLExporter>::New();
  exporter->SetRenderWindow(_rw1);
  exporter->SetFileName((path+".vrml").c_str());
  exporter->Write();
}

void G4VtkViewer::ExportVTPScene(G4String path)
{
  vtkSmartPointer<vtkRenderWindow> _rw1 = vtkSmartPointer<vtkRenderWindow>::New();
  _rw1->AddRenderer(_renderWindow->GetRenderers()->GetFirstRenderer());
  vtkSmartPointer<vtkSingleVTPExporter> exporter = vtkSmartPointer<vtkSingleVTPExporter>::New();
  exporter->SetRenderWindow(_rw1);
  exporter->SetFileName((path+".vtp").c_str());
  exporter->Write();
}
