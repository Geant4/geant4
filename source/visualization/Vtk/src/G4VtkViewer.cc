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

#include <cmath>

#include "G4VtkViewer.hh"

#include "G4Transform3D.hh"
#include "G4VSceneHandler.hh"
#include "G4VtkCutterPipeline.hh"
#include "G4VtkPolydataInstanceAppendPipeline.hh"
#include "G4VtkPolydataInstanceBakePipeline.hh"
#include "G4VtkPolydataInstancePipeline.hh"
#include "G4VtkPolydataInstanceTensorPipeline.hh"
#include "G4VtkSceneHandler.hh"
#include "G4VtkUtility.hh"
#include "G4VtkVisContext.hh"

#include "vtk3DSImporter.h"
#include "vtkBMPWriter.h"
#include "vtkIVExporter.h"  // open inventor
#include "vtkImageWriter.h"
#include "vtkImplicitPlaneRepresentation.h"
#include "vtkImplicitPlaneWidget2.h"
#include "vtkJPEGWriter.h"
#include "vtkLightCollection.h"
#include "vtkOBJExporter.h"
#include "vtkOBJImporter.h"
#include "vtkGLTFExporter.h"
#include "vtkOOGLExporter.h"
#include "vtkJSONRenderWindowExporter.h"
#include "vtkVtkJSSceneGraphSerializer.h"
//#include "vtkBufferedArchiver.h"
#include "vtkPNGWriter.h"
#include "vtkPNMWriter.h"
#include "vtkPOVExporter.h"
#include "vtkPostScriptWriter.h"
#include "vtkRIBExporter.h"  // Renderman
#include "vtkRendererCollection.h"
#include "vtkSingleVTPExporter.h"
#include "vtkTIFFWriter.h"
#include "vtkVRMLExporter.h"
#include "vtkVRMLImporter.h"
#include "vtkWindowToImageFilter.h"
#include "vtkX3DExporter.h"

// Readers (vtkDataReader)

#include "vtkCameraPass.h"
#include "vtkOpenGLRenderer.h"
#include "vtkRenderPass.h"
#include "vtkRenderPassCollection.h"
#include "vtkSequencePass.h"
#include "vtkShadowMapBakerPass.h"
#include "vtkShadowMapPass.h"

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
  _renderWindow = vtkRenderWindow::New();
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
  _renderWindow->SetSize(pixelFactor * fVP.GetWindowSizeHintX(),
                         pixelFactor * fVP.GetWindowSizeHintY());
  _renderWindow->SetWindowName("Vtk viewer");

  _renderWindow->AddRenderer(renderer);
  renderWindowInteractor->SetRenderWindow(_renderWindow);

  // TODO proper camera parameter settings
  camera->SetPosition(0, 0, 1000);
  camera->SetFocalPoint(0, 0, 0);
  renderer->SetActiveCamera(camera);

  // Hidden line removal
  renderer->SetUseHiddenLineRemoval(0);

  // Shadows
  renderer->SetUseShadows(0);

  // Set callback to match VTK parameters to Geant4
  geant4Callback->SetGeant4ViewParameters(&fVP);
  renderer->AddObserver(vtkCommand::EndEvent, geant4Callback);

  vtkSmartPointer<G4VtkInteractorStyle> style = vtkSmartPointer<G4VtkInteractorStyle>::New();
  renderWindowInteractor->SetInteractorStyle(style);
}

G4VtkViewer::~G4VtkViewer()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::~G4VtkViewer()" << G4endl;
#endif
}

void G4VtkViewer::SetView()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::SetView()" << G4endl;
#endif
  // background colour
  const G4Colour backgroundColour = fVP.GetBackgroundColour();
  renderer->SetBackground(backgroundColour.GetRed(), backgroundColour.GetGreen(),
                          backgroundColour.GetBlue());

  // target and camera positions
  G4double radius = fSceneHandler.GetExtent().GetExtentRadius();
  if (radius <= 0.) {
    radius = 1.;
  }

  if (firstSetView) cameraDistance = fVP.GetCameraDistance(radius);
  G4Point3D viewpointDirection = fVP.GetViewpointDirection();
  G4Point3D targetPosition = fVP.GetCurrentTargetPoint();
  G4double zoomFactor = fVP.GetZoomFactor();
  G4double fieldHalfAngle = fVP.GetFieldHalfAngle();

  G4Point3D cameraPosition = targetPosition + viewpointDirection.unit() /zoomFactor  * cameraDistance;

  vtkCamera* activeCamera = renderer->GetActiveCamera();
  activeCamera->SetFocalPoint(targetPosition.x(), targetPosition.y(), targetPosition.z());
  activeCamera->SetViewAngle(2*fieldHalfAngle / M_PI * 180);
  activeCamera->SetPosition(cameraPosition.x(), cameraPosition.y(), cameraPosition.z());
  activeCamera->SetParallelScale(cameraDistance / zoomFactor);

  if (fieldHalfAngle == 0) {
    activeCamera->SetParallelProjection(1);
  }
  else {
    activeCamera->SetParallelProjection(0);
  }

  // need to set camera distance and parallel scale on first set view
  if (firstSetView) {
    geant4Callback->SetVtkInitialValues(cameraDistance, cameraDistance);
    activeCamera->SetParallelScale(cameraDistance);
    firstSetView = false;
  }

  // camera up direction
  const G4Vector3D upVector = fVP.GetUpVector();
  renderer->GetActiveCamera()->SetViewUp(upVector.x(), upVector.y(), upVector.z());

  // Light
  const G4Vector3D lightDirection = fVP.GetLightpointDirection();
  G4bool lightsMoveWithCamera = fVP.GetLightsMoveWithCamera();
  G4Vector3D lightPosition = targetPosition + lightDirection.unit() * cameraDistance;

  vtkLightCollection* currentLights = renderer->GetLights();
  if (currentLights->GetNumberOfItems() != 0) {
    auto currentLight = dynamic_cast<vtkLight*>(currentLights->GetItemAsObject(0));
    if (currentLight != nullptr) {
      currentLight->SetPosition(lightPosition.x(), lightPosition.y(), lightPosition.z());
      if (lightsMoveWithCamera) {
        currentLight->SetLightTypeToCameraLight();
      }
      else {
        currentLight->SetLightTypeToSceneLight();
      }
    }
  }

  // cut away
  if (fVP.IsCutaway()) {
    G4cout << "Add cutaway planes" << G4endl;
  }

  // section
  if (fVP.IsSection()) {
    G4cout << "Add section" << G4endl;
  }
}

void G4VtkViewer::ClearView()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ClearView()" << G4endl;
#endif

  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& ts = fVtkSceneHandler.GetTransientStore();
  ts.Clear();

  G4VtkStore& s = fVtkSceneHandler.GetStore();
  s.Clear();
}

void G4VtkViewer::DrawView()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::DrawView()" << G4endl;
#endif

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
  AddViewHUD();

  // Add clipper and cutter widgets
  auto g4p = G4Plane3D();
  AddCutterPlaneWidget(g4p);
  AddClipperPlaneWidget(g4p);

  // Add camera orientation widget
  AddCameraOrientationWidget();

  // ...before finally...
  FinishView();  // Flush streams and/or swap buffers.
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
  auto glrenderer = dynamic_cast<vtkOpenGLRenderer*>(renderer.GetPointer());
  glrenderer->SetPass(cameraP);
}

void G4VtkViewer::ShowView()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ShowView()" << G4endl;
#endif

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
}

void G4VtkViewer::FinishView()
{
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::FinishView()" << G4endl;
#endif

  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Modified();

  _renderWindow->GetInteractor()->Initialize();
  _renderWindow->Render();

  if (firstFinishView) {
    firstFinishView = false;
  }
  else {
    _renderWindow->GetInteractor()->Start();
  }
}

void G4VtkViewer::ExportScreenShot(G4String path, G4String format)
{
  vtkImageWriter* imWriter = nullptr;

  if (format == "bmp") {
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
    return;
  }

  _renderWindow->Render();

  vtkSmartPointer<vtkWindowToImageFilter> winToImage =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  winToImage->SetInput(_renderWindow);
  winToImage->SetScale(1);
  if (format == "ps") {
    winToImage->SetInputBufferTypeToRGB();
    winToImage->ReadFrontBufferOff();
    winToImage->Update();
  }
  else {
    winToImage->SetInputBufferTypeToRGBA();
  }

  imWriter->SetFileName((path + "." + format).c_str());
  imWriter->SetInputConnection(winToImage->GetOutputPort());
  imWriter->Write();
}

void G4VtkViewer::ExportOBJScene(G4String path)
{
  vtkSmartPointer<vtkOBJExporter> exporter = vtkSmartPointer<vtkOBJExporter>::New();
  exporter->SetRenderWindow(_renderWindow);
  exporter->SetFilePrefix(path.c_str());
  exporter->Write();
}

void G4VtkViewer::ExportVRMLScene(G4String path)
{
  vtkSmartPointer<vtkVRMLExporter> exporter = vtkSmartPointer<vtkVRMLExporter>::New();
  exporter->SetRenderWindow(_renderWindow);
  exporter->SetFileName((path + ".vrml").c_str());
  exporter->Write();
}

void G4VtkViewer::ExportVTPScene(G4String path)
{
  vtkSmartPointer<vtkSingleVTPExporter> exporter = vtkSmartPointer<vtkSingleVTPExporter>::New();
  exporter->SetRenderWindow(_renderWindow);
  exporter->SetFileName((path + ".vtp").c_str());
  exporter->Write();
}

void G4VtkViewer::ExportGLTFScene(G4String fileName) {
  vtkSmartPointer<vtkGLTFExporter> exporter = vtkSmartPointer<vtkGLTFExporter>::New();
  exporter->SetRenderWindow(_renderWindow);
  exporter->SetFileName((fileName+".gltf").c_str());
  exporter->InlineDataOn();
  exporter->Write();
}

void G4VtkViewer::ExportJSONRenderWindowScene(G4String /*fileName*/) {
  vtkSmartPointer<vtkJSONRenderWindowExporter> exporter = vtkSmartPointer<vtkJSONRenderWindowExporter>::New();
  vtkSmartPointer<vtkVtkJSSceneGraphSerializer> serializer = vtkSmartPointer<vtkVtkJSSceneGraphSerializer>::New();
  exporter->SetRenderWindow(_renderWindow);
  exporter->SetSerializer(serializer);
  exporter->Write();
}

void G4VtkViewer::ExportVTPCutter(G4String fileName)
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& s = fVtkSceneHandler.GetStore();

  // create new renderer
  vtkNew<vtkRenderer> tempRenderer;

  // loop over pipelines
  auto separate = s.GetSeparatePipeMap();
  for (const auto& i : separate) {
    i.second->GetActor();
    auto children = i.second->GetChildPipelines();
    for (auto child : children) {
      if (child->GetTypeName() == "G4VtkCutterPipeline") {
        auto childCutter = dynamic_cast<G4VtkCutterPipeline*>(child);
        tempRenderer->AddActor(childCutter->GetActor());
      }
    }
  }

  auto tensor = s.GetTensorPipeMap();
  for (const auto& i : tensor) {
    i.second->GetActor();
    auto children = i.second->GetChildPipelines();
    for (auto child : children) {
      if (child->GetTypeName() == "G4VtkCutterPipeline") {
        auto childCutter = dynamic_cast<G4VtkCutterPipeline*>(child);
        tempRenderer->AddActor(childCutter->GetActor());
      }
    }
  }

  auto append = s.GetAppendPipeMap();
  for (const auto& i : append) {
    i.second->GetActor();
    auto children = i.second->GetChildPipelines();
    for (auto child : children) {
      if (child->GetTypeName() == "G4VtkCutterPipeline") {
        auto childCutter = dynamic_cast<G4VtkCutterPipeline*>(child);
        tempRenderer->AddActor(childCutter->GetActor());
      }
    }
  }

  auto baked = s.GetBakePipeMap();
  for (const auto& i : baked) {
    i.second->GetActor();
    auto children = i.second->GetChildPipelines();
    for (auto child : children) {
      if (child->GetTypeName() == "G4VtkCutterPipeline") {
        auto childCutter = dynamic_cast<G4VtkCutterPipeline*>(child);
        tempRenderer->AddActor(childCutter->GetActor());
      }
    }
  }

  vtkNew<vtkRenderWindow> tempRenderWindow;
  tempRenderWindow->AddRenderer(tempRenderer);
  vtkNew<vtkSingleVTPExporter> exporter;
  exporter->SetRenderWindow(tempRenderWindow);
  exporter->SetFileName(fileName.c_str());
  exporter->Write();
}

void G4VtkViewer::ExportFormatStore(G4String fileName, G4String storeName)
{
  vtkSmartPointer<vtkRenderWindow> tempRenderWindow;
  vtkNew<vtkRenderer> tempRenderer;
  tempRenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  tempRenderWindow->AddRenderer(tempRenderer);

  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);

  if (storeName == "transient") {
    G4VtkStore& store = fVtkSceneHandler.GetTransientStore();
    store.AddToRenderer(tempRenderer);
  }
  else {
    G4VtkStore& store = fVtkSceneHandler.GetStore();
    store.AddToRenderer(tempRenderer);
  }

  if (fileName.find("obj") != std::string::npos) {
    vtkNew<vtkOBJExporter> exporter;
    exporter->SetRenderWindow(tempRenderWindow);
    exporter->SetFilePrefix(fileName.c_str());
    exporter->Write();
  }
  else if (fileName.find("vrml") != std::string::npos) {
    vtkNew<vtkVRMLExporter> exporter;
    exporter->SetRenderWindow(tempRenderWindow);
    exporter->SetFileName(fileName.c_str());
    exporter->Write();
  }
  else if (fileName.find("vtp") != std::string::npos) {
    vtkNew<vtkSingleVTPExporter> exporter;
    exporter->SetRenderWindow(tempRenderWindow);
    exporter->SetFileName(fileName.c_str());
    exporter->Write();
  }
  else if (fileName.find("gltf") != std::string::npos) {
    vtkNew<vtkGLTFExporter> exporter;
    exporter->SetRenderWindow(tempRenderWindow);
    exporter->SetFileName(fileName.c_str());
    exporter->InlineDataOn();
    exporter->Write();
  }
}

void G4VtkViewer::AddViewHUD()
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
  infoTextActor->SetVisibility(0);
}

void G4VtkViewer::AddClipperPlaneWidget(const G4Plane3D& plane)
{
  vtkNew<vtkIPWCallback> clipperCallback;
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& store = fVtkSceneHandler.GetStore();
  clipperCallback->SetStore(&store);
  clipperCallback->SetUpdatePipelineName("clipper", "clipper");

  G4double bounds[6];
  store.GetBounds(bounds);
  auto vplane = G4Plane3DToVtkPlane(plane);
  clipperPlaneRepresentation->SetPlaceFactor(
    1.25);  // This must be set prior to placing the widget.
  clipperPlaneRepresentation->PlaceWidget(bounds);
  clipperPlaneRepresentation->SetNormal(vplane->GetNormal());

  vtkNew<vtkPropCollection> planeRepActors;
  clipperPlaneRepresentation->GetActors(planeRepActors);
  planeRepActors->InitTraversal();

  SetWidgetInteractor(clipperPlaneWidget);
  clipperPlaneWidget->SetRepresentation(clipperPlaneRepresentation);
  clipperPlaneWidget->AddObserver(vtkCommand::InteractionEvent, clipperCallback);

  clipperPlaneWidget->SetEnabled(0);
}

void G4VtkViewer::AddCutterPlaneWidget(const G4Plane3D& plane)
{
  vtkNew<vtkIPWCallback> cutterCallback;
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& store = fVtkSceneHandler.GetStore();
  cutterCallback->SetStore(&store);
  cutterCallback->SetUpdatePipelineName("cutter", "cutter");

  G4double bounds[6];
  store.GetBounds(bounds);
  auto vplane = G4Plane3DToVtkPlane(plane);
  cutterPlaneRepresentation->SetPlaceFactor(1.25);  // This must be set prior to placing the widget.
  cutterPlaneRepresentation->PlaceWidget(bounds);
  cutterPlaneRepresentation->SetNormal(vplane->GetNormal());

  SetWidgetInteractor(cutterPlaneWidget);
  cutterPlaneWidget->SetRepresentation(cutterPlaneRepresentation);
  cutterPlaneWidget->AddObserver(vtkCommand::InteractionEvent, cutterCallback);

  cutterPlaneWidget->SetEnabled(0);
}

void G4VtkViewer::EnableShadows()
{
  renderer->SetUseShadows(1);
}

void G4VtkViewer::DisableShadows()
{
  renderer->SetUseShadows(0);
}

void G4VtkViewer::EnableHUD()
{
  infoTextActor->SetVisibility(1);
}

void G4VtkViewer::DisableHUD()
{
  infoTextActor->SetVisibility(0);
}

void G4VtkViewer::EnableClipper(const G4Plane3D& plane, G4bool bWidget)
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& s = fVtkSceneHandler.GetStore();
  G4String name = G4String("clipper");
  s.AddClipper(name, plane);
  if (bWidget) {
    EnableClipperWidget();
  }
}

void G4VtkViewer::DisableClipper()
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& s = fVtkSceneHandler.GetStore();
  s.RemoveClipper("clipper");
}

void G4VtkViewer::EnableClipperWidget()
{
  clipperPlaneWidget->SetEnabled(1);
}

void G4VtkViewer::DisableClipperWidget()
{
  clipperPlaneWidget->SetEnabled(0);
}

void G4VtkViewer::EnableCutter(const G4Plane3D& plane, G4bool bWidget)
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& s = fVtkSceneHandler.GetStore();
  G4String name = G4String("cutter");
  s.AddCutter(name, plane);
  if (bWidget) {
    EnableCutterWidget();
  }
}

void G4VtkViewer::DisableCutter(G4String /*name*/)
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& s = fVtkSceneHandler.GetStore();
  s.RemoveCutter("cutter");
}

void G4VtkViewer::EnableCutterWidget()
{
  G4cout << "enable cutter widget" << G4endl;
  cutterPlaneWidget->SetEnabled(1);
}

void G4VtkViewer::DisableCutterWidget()
{
  cutterPlaneWidget->SetEnabled(0);
}

void G4VtkViewer::AddCameraOrientationWidget()
{
  camOrientWidget->SetParentRenderer(renderer);
  // Enable the widget.
  camOrientWidget->Off();
}

void G4VtkViewer::EnableCameraOrientationWidget()
{
  camOrientWidget->On();
}

void G4VtkViewer::DisableCameraOrientationWidget()
{
  camOrientWidget->Off();
}

void G4VtkViewer::AddImageOverlay(const G4String& fileName, const G4double alpha,
                                  const G4double imageBottomLeft[2],
                                  const G4double worldBottomLeft[2],
                                  const G4double imageTopRight[2], const G4double worldTopRight[2],
                                  const G4double rotation[3], const G4double translation[3])
{
  auto xScale = (worldTopRight[0] - worldBottomLeft[0]) / (imageTopRight[0] - imageBottomLeft[0]);
  auto yScale = -(worldTopRight[1] - worldBottomLeft[1]) / (imageTopRight[1] - imageBottomLeft[1]);

  G4cout << xScale << " " << yScale << G4endl;
  auto transformation = G4Transform3D::Identity;
  auto scal = G4Scale3D(xScale, yScale, 1);
  auto rotx = G4RotateX3D(rotation[0]/180*M_PI);
  auto roty = G4RotateY3D(rotation[1]/180*M_PI);
  auto rotz = G4RotateZ3D(rotation[2]/180*M_PI);
  auto tranImg = G4Translate3D(  -std::fabs(imageBottomLeft[0] + imageTopRight[0]) / 2.0,
                                 -std::fabs(imageBottomLeft[1] + imageTopRight[1]) / 2.0,
                                 0);
  auto tran = G4Translate3D(translation[0],
                            translation[1],
                            translation[2]);

  G4cout << translation[0] << " " << translation[1] << " " << translation[2] << G4endl;
  transformation = tran * rotz * roty * rotx * scal * tranImg * transformation;

  G4cout << transformation.dx() << " " << transformation.dy() << " " << transformation.dz() << G4endl;
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& st = fVtkSceneHandler.GetTransientStore();

  G4VtkVisContext vc = G4VtkVisContext(this, nullptr, false, transformation);
  vc.alpha = alpha;
  st.AddNonG4ObjectImage(fileName, vc);
}

void G4VtkViewer::AddGeometryOverlay(const G4String& fileName, const G4double colour[3],
                                     const G4double alpha, const G4String& representation,
                                     const G4double scale[3], const G4double rotation[3],
                                     const G4double translation[3])
{
  auto transformation = G4Transform3D::Identity;
  auto scal = G4Scale3D(scale[0], scale[1], scale[2]);
  auto rotx = G4RotateX3D(rotation[0]);
  auto roty = G4RotateY3D(rotation[1]);
  auto rotz = G4RotateZ3D(rotation[2]);
  auto tran = G4Translate3D(translation[0], translation[1], translation[2]);

  transformation = tran * rotz * roty * rotx * scal * transformation;

  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& st = fVtkSceneHandler.GetTransientStore();


  G4VtkVisContext vc = G4VtkVisContext(this, nullptr, false, transformation);
  if (representation == "w")
    vc.fDrawingStyle = G4ViewParameters::wireframe;
  else if (representation == "s")
    vc.fDrawingStyle = G4ViewParameters::hlhsr;
  vc.alpha = alpha;
  vc.red = colour[0];
  vc.green = colour[1];
  vc.blue = colour[2];
  st.AddNonG4ObjectPolydata(fileName, vc);
}

void G4VtkViewer::Print()
{
  cutterPlaneRepresentation->VisibilityOff();

  G4cout << "Number of VTK props>  " << renderer->GetNumberOfPropsRendered() << G4endl;
  G4cout << "Number of VTK actors> " << renderer->GetActors()->GetNumberOfItems() << G4endl;
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  G4VtkStore& s = fVtkSceneHandler.GetStore();
  G4VtkStore& st = fVtkSceneHandler.GetTransientStore();
  s.Print();
  st.Print();
}

void G4VtkViewer::SetPolyhedronPipeline(const G4String& type)
{
  // Get the scene handler
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.SetPolyhedronPipeline(type);
}

void G4VtkViewer::SetWidgetInteractor(vtkAbstractWidget* widget)
{
  widget->SetInteractor(_renderWindow->GetInteractor());
}
