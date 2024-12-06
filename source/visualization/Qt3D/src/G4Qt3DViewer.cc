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
// John Allison  17th June 2019

#include "G4Qt3DViewer.hh"

#include "G4Qt3DSceneHandler.hh"
#include "G4Qt3DUtils.hh"

#include "G4Scene.hh"
#include "G4UImanager.hh"
#include "G4UIQt.hh"
#include "G4SystemOfUnits.hh"

#define G4warn G4cout

G4Qt3DViewer::G4Qt3DViewer
(G4Qt3DSceneHandler& sceneHandler, const G4String& name)
: G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
, fQt3DSceneHandler(sceneHandler)
, fKeyPressed(false)
, fMousePressed(false)
, fMousePressedX(0.)
, fMousePressedY(0.)
{}

void G4Qt3DViewer::Initialise()
{
  setObjectName(fName.c_str());

  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);

  auto UI = G4UImanager::GetUIpointer();
  auto uiQt = dynamic_cast<G4UIQt*>(UI->GetG4UIWindow());
  if (!uiQt) {
    fViewId = -1;  // This flags an error.
    G4warn << "G4Qt3DViewer::G4Qt3DViewer requires G4UIQt"
         << G4endl;
    return;
  }
  fUIWidget = QWidget::createWindowContainer(this);
  uiQt->AddTabWidget(fUIWidget,QString(fName));

  setRootEntity(fQt3DSceneHandler.fpQt3DScene);
}

G4Qt3DViewer::~G4Qt3DViewer()
{
  setRootEntity(nullptr);
}

void G4Qt3DViewer::resizeEvent(QResizeEvent*) {
  SetView();
}

void G4Qt3DViewer::SetView()
{
  // Background colour
  defaultFrameGraph()->setClearColor(G4Qt3DUtils::ConvertToQColor(fVP.GetBackgroundColour()));

  // Get radius of scene, etc.
  // Note that this procedure properly takes into account zoom, dolly and pan.
  const G4Point3D targetPoint
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const G4Point3D cameraPosition =
    targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
  const GLdouble pnear  = fVP.GetNearDistance (cameraDistance, radius);
  const GLdouble pfar   = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  const GLdouble right  = fVP.GetFrontHalfHeight (pnear, radius);
  const GLdouble left   = -right;
  const GLdouble top    = fVP.GetFrontHalfHeight (pnear, radius);
  const GLdouble bottom = -top;

  camera()->setObjectName((fName + " camera").c_str());
  camera()->setViewCenter(G4Qt3DUtils::ConvertToQVector3D(targetPoint));
  camera()->setPosition(G4Qt3DUtils::ConvertToQVector3D(cameraPosition));
  camera()->setUpVector(G4Qt3DUtils::ConvertToQVector3D(fVP.GetUpVector()));

//  auto lightEntity = new Qt3DCore::QEntity(fQt3DSceneHandler.fpQt3DScene);
//  auto directionalLight = new Qt3DRender::QDirectionalLight(lightEntity);
////  directionalLight->setColor("white");
////  directionalLight->setIntensity(1.);
//  directionalLight->setWorldDirection(G4Qt3DUtils::ConvertToQVector3D(fVP.GetActualLightpointDirection()));
//  lightEntity->addComponent(directionalLight);

  const auto& size = fUIWidget->size();
  G4double w = size.width();
  G4double h = size.height();
#ifdef G4QT3DDEBUG
  // Curiously w,h are wrong first time - 640,480 instead of (my Mac) 991,452.
  G4cout << "W,H: " << w << ',' << h << G4endl;
#endif
  const G4double aspectRatio = w/h;
  if (fVP.GetFieldHalfAngle() == 0.) {
    camera()->lens()->setOrthographicProjection
    (left*aspectRatio,right*aspectRatio,bottom,top,pnear,pfar);
  } else {
    camera()->lens()->setPerspectiveProjection
    (2.*fVP.GetFieldHalfAngle()/deg,aspectRatio,pnear,pfar);
  }
}

void G4Qt3DViewer::ClearView()
{}

void G4Qt3DViewer::DrawView()
{
  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  // The fNeedKernelVisit flag might have been set by the user in
  // /vis/viewer/rebuild, but if not, make decision and set flag only
  // if necessary...
  if (!fNeedKernelVisit) KernelVisitDecision();
  G4bool kernelVisitWasNeeded = fNeedKernelVisit; // Keep (ProcessView resets).
  fLastVP = fVP;

  ProcessView ();  // Clears store and processes scene only if necessary.

  if (kernelVisitWasNeeded) {
    // We might need to do something if the kernel was visited.
  } else {
  }

  // ...before finally...
  FinishView ();       // Flush streams and/or swap buffers.
}

void G4Qt3DViewer::ShowView()
{
#if QT_VERSION < 0x060000
  // show() may only be called from master thread
  if (G4Threading::IsMasterThread()) {
    show();
  }
  // The way Qt seems to work, we don't seem to need a show() anyway, but
  // we'll leave it in - it seems not to have any effect, good or bad.
#endif
}

void G4Qt3DViewer::FinishView()
{
#if QT_VERSION < 0x060000
  if (G4Threading::IsMasterThread()) {
    show();
  }
#endif
}

// Note: the order of calling of MovingToVisSubThread and SwitchToVisSubThread
// is undefined. The order of calling is
//   DoneWithMasterThread
//   MovingToVisSubThread ) or ( SwitchToVisSubThread
//   SwitchToVisSubThread )    ( MovingToVisSubThread
//   DoneWithVisSubThread
//   MovingToMasterThread
//   SwitchToMasterThread
// So regarding the move/switch to the vis sub-thread, we have to employ mutex locks and conditions.
// If the viewer wishes to accept drawing from the vis sub-thread, Qt3D has to moveToThread.
// But at this point we are still on the master thread and the value of the sub-thread's QThread is
// not known. So it has to wait - a conditional wait, conditional on establishment of the sub-thread
// and the provision of a pointer to the QThread version of the vis sub-thread. In turn, the
// sub-thread has to wait until the QObjects have been moved to the sub-thread.

namespace {
  QThread* masterQThread = nullptr;
  QThread* visSubThreadQThread = nullptr;

  G4Mutex visSubThreadMutex = G4MUTEX_INITIALIZER;
  G4Condition waitForVisSubThreadInitialized = G4CONDITION_INITIALIZER;
  G4bool visSubThreadEstablished = false;
  G4bool qObjectsSwitched = false;
}

void G4Qt3DViewer::MovingToVisSubThread()
// Still on master thread but vis thread has been launched
{
  // Make note of master QThread
  masterQThread = QThread::currentThread();

  // Wait until SwitchToVisSubThread has found vis sub-thread QThread
  {
  G4AutoLock lock(&visSubThreadMutex);
  G4CONDITIONWAITLAMBDA(&waitForVisSubThreadInitialized, &lock, []{return visSubThreadEstablished;})
  }

  // Move relevant stuff to vis sub-thread QThread
  auto p1 = fQt3DSceneHandler.fpQt3DScene->parent();
  if(p1) {
    auto p2 = p1->parent();
    if(p2) {
      p2->moveToThread(visSubThreadQThread);
    } else {
      p1->moveToThread(visSubThreadQThread);
    }
  }

  // Inform sub-thread
  G4AutoLock lock(&visSubThreadMutex);
  qObjectsSwitched = true;
  lock.unlock();
  G4CONDITIONBROADCAST(&waitForVisSubThreadInitialized);
}

void G4Qt3DViewer::SwitchToVisSubThread()
// On vis sub-thread before any drawing
{
  // Make note of vis-subthread QThread for MovingToVisSubThread
  visSubThreadQThread = QThread::currentThread();

  // Let MovingToVisSubThread know we have the QThread
  {
  G4AutoLock lock(&visSubThreadMutex);
  visSubThreadEstablished = true;
  G4CONDITIONBROADCAST(&waitForVisSubThreadInitialized);
  }

  // Wait until MovingToVisSubThread has moved stuff
  {
  G4AutoLock lock(&visSubThreadMutex);
  G4CONDITIONWAITLAMBDA(&waitForVisSubThreadInitialized, &lock, []{return qObjectsSwitched;})
  }
}

void G4Qt3DViewer::MovingToMasterThread()
// On vis sub-thread just before exit
{
  // Move relevant stuff to master QThread.
  auto p1 = fQt3DSceneHandler.fpQt3DScene->parent();
  if(p1) {
    auto p2 = p1->parent();
    if(p2) {
      p2->moveToThread(masterQThread);
    } else {
      p1->moveToThread(masterQThread);
    }
  }

  // Reset
  visSubThreadQThread = nullptr;
  qObjectsSwitched = false;
}

void G4Qt3DViewer::SwitchToMasterThread()
// On master thread after vis sub-thread has terminated
{
  visSubThreadEstablished = false;
}

void G4Qt3DViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  if (CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();  // Sets fNeedKernelVisit.
  }
}

G4bool G4Qt3DViewer::CompareForKernelVisit(G4ViewParameters& vp)
{
  // Typical comparison.  Taken from OpenInventor.
  if (
     (vp.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
     (vp.GetNumberOfCloudPoints()  != fVP.GetNumberOfCloudPoints())  ||
     (vp.IsAuxEdgeVisible ()   != fVP.IsAuxEdgeVisible ())   ||
     (vp.IsCulling ()          != fVP.IsCulling ())          ||
     (vp.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
     (vp.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
     (vp.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
     (vp.GetCBDAlgorithmNumber() !=
      fVP.GetCBDAlgorithmNumber())                           ||
     (vp.IsSection ()          != fVP.IsSection ())          ||
     (vp.IsCutaway ()          != fVP.IsCutaway ())          ||
     // This assumes use of generic clipping (sectioning, slicing,
     // DCUT, cutaway).  If a decision is made to implement locally,
     // this will need changing.  See G4OpenGLViewer::SetView,
     // G4OpenGLStoredViewer.cc::CompareForKernelVisit and
     // G4OpenGLStoredSceneHander::CreateSection/CutawayPolyhedron.
     (vp.IsExplode ()          != fVP.IsExplode ())          ||
     (vp.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
     (vp.GetGlobalMarkerScale()    != fVP.GetGlobalMarkerScale())    ||
     (vp.GetGlobalLineWidthScale() != fVP.GetGlobalLineWidthScale()) ||
     (vp.IsMarkerNotHidden ()  != fVP.IsMarkerNotHidden ())  ||
     (vp.GetDefaultVisAttributes()->GetColour() !=
      fVP.GetDefaultVisAttributes()->GetColour())            ||
     (vp.GetDefaultTextVisAttributes()->GetColour() !=
      fVP.GetDefaultTextVisAttributes()->GetColour())        ||
     (vp.GetBackgroundColour ()!= fVP.GetBackgroundColour ())||
     (vp.IsPicking ()          != fVP.IsPicking ())          ||
     // Scaling for Open Inventor is done by the scene handler so it
     // needs a kernel visit.  (In this respect, it differs from the
     // OpenGL drivers, where it's done in SetView.)
     (vp.GetScaleFactor ()     != fVP.GetScaleFactor ())     ||
     (vp.GetVisAttributesModifiers() !=
      fVP.GetVisAttributesModifiers())                       ||
     (vp.IsSpecialMeshRendering() !=
      fVP.IsSpecialMeshRendering())                          ||
     (vp.GetSpecialMeshRenderingOption() !=
      fVP.GetSpecialMeshRenderingOption())
     )
  return true;

  if (vp.IsDensityCulling () &&
      (vp.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  if (vp.GetCBDAlgorithmNumber() > 0) {
    if (vp.GetCBDParameters().size() != fVP.GetCBDParameters().size()) return true;
    else if (vp.GetCBDParameters() != fVP.GetCBDParameters()) return true;
  }

  if (vp.IsSection () &&
      (vp.GetSectionPlane () != fVP.GetSectionPlane ()))
    return true;

  if (vp.IsCutaway ()) {
    if (vp.GetCutawayMode() != fVP.GetCutawayMode()) return true;
    if (vp.GetCutawayPlanes ().size () !=
        fVP.GetCutawayPlanes ().size ()) return true;
    for (size_t i = 0; i < vp.GetCutawayPlanes().size(); ++i)
    if (vp.GetCutawayPlanes()[i] != fVP.GetCutawayPlanes()[i])
      return true;
  }

  if (vp.IsExplode () &&
      (vp.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;

  if (vp.IsSpecialMeshRendering() &&
      (vp.GetSpecialMeshVolumes() != fVP.GetSpecialMeshVolumes()))
    return true;

  return false;
}

void G4Qt3DViewer::keyPressEvent(QKeyEvent* ev)
{
  fKeyPressed = true;
  fKey = ev->key();
}

void G4Qt3DViewer::keyReleaseEvent(QKeyEvent* /*ev*/)
{
  fKeyPressed = false;
}

void G4Qt3DViewer::mouseDoubleClickEvent(QMouseEvent* /*ev*/) {}

void G4Qt3DViewer::mouseMoveEvent(QMouseEvent* ev)
{
  // I think we only want these if a mouse button is pressed.
  // But they come even when not pressed (on my MacBook Pro trackpad).
  // Documentation says:
  /* Mouse move events will occur only when a mouse button is pressed down,
   unless mouse tracking has been enabled with QWidget::setMouseTracking().*/
  // But this is a window not a widget.
  // As a workaround we maintain a flag changed by mousePress/ReleaseEvent.
#if (QT_VERSION < QT_VERSION_CHECK(6, 0, 0))
  G4double x = ev->x();
  G4double y = ev->y();
#else
  G4double x = ev->position().x();
  G4double y = ev->position().y();
#endif
  G4double dx = x-fMousePressedX;
  G4double dy = y-fMousePressedY;
  fMousePressedX = x;
  fMousePressedY = y;

  if (fMousePressed) {

    if (fKeyPressed && fKey == Qt::Key_Shift) {  // Translation (pan)

      const G4double sceneRadius = fQt3DSceneHandler.fpScene->GetExtent().GetExtentRadius();
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

void G4Qt3DViewer::mousePressEvent(QMouseEvent* ev)
{
  fMousePressed = true;
#if (QT_VERSION < QT_VERSION_CHECK(6, 0, 0))
  fMousePressedX = ev->x();
  fMousePressedY = ev->y();
#else
  fMousePressedX = ev->position().x();
  fMousePressedY = ev->position().y();
#endif
}

void G4Qt3DViewer::mouseReleaseEvent(QMouseEvent* /*ev*/)
{
  fMousePressed = false;
}

void G4Qt3DViewer::wheelEvent(QWheelEvent* ev)
{
  // Take note of up-down motion only
  const G4double angleY = ev->angleDelta().y();

  if (fVP.GetFieldHalfAngle() == 0.) {  // Orthographic projection
    const G4double scale = 500;  // Empirically chosen
    fVP.MultiplyZoomFactor(1.+angleY/scale);
  } else {                              // Perspective projection
    const G4double delta = fSceneHandler.GetExtent().GetExtentRadius()/200.;  // Empirical
    fVP.SetDolly(fVP.GetDolly()+angleY*delta);
  }
  
  SetView();
  DrawView();
}
