//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "G4OpenInventorViewer.hh"

#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoShape.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/sensors/SoNodeSensor.h>

#include "HEPVis/nodes/SoImageWriter.h"
#include "HEPVis/actions/SoGL2PSAction.h"
#include "HEPVis/actions/SoCounterAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"
#include "G4Scene.hh"
#include "Geant4_SoPolyhedron.h"

G4OpenInventorViewer::G4OpenInventorViewer(
 G4OpenInventorSceneHandler& sceneHandler
,const G4String& name)
:G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
,fG4OpenInventorSceneHandler(sceneHandler)
,fInteractorManager(0)
,fSoSelection(0)
,fSoCamera(0)
,fSoImageWriter(0)
,fGL2PSAction(0) //To be set be suclass.
,fSoCameraSensor(0)
{
  fNeedKernelVisit = true;  //?? Temporary, until KernelVisitDecision fixed.

  fDefaultVP.SetAutoRefresh(true);
  fVP.SetAutoRefresh(true);

  //FIXME : G.Barrand : not convinced that we have to rm culling.
  // For viewing of all objects by default :
  //fDefaultVP.SetCulling(false);
  //fVP.SetCulling(false);

  fInteractorManager = 
    ((G4OpenInventor*)fG4OpenInventorSceneHandler.GetGraphicsSystem())->
    GetInteractorManager();

  // Main user scene graph root sent to the viewers.
  fSoSelection = new SoSelection;

  fSoSelection->addSelectionCallback(SelectionCB,this);
  //fSoSelection->addDeselectionCallback(DeselectionCB,this);

  fSoSelection->policy = SoSelection::SINGLE;
  fSoSelection->ref();

  fSoCamera = new SoOrthographicCamera;
  fSoCamera->viewportMapping.setValue(SoCamera::ADJUST_CAMERA);
  //camera->aspectRatio.setValue(10);
  fSoCamera->position.setValue(0,0,10);
  fSoCamera->orientation.setValue(SbRotation(SbVec3f(0,1,0),0));
  fSoCamera->height.setValue(10);
  fSoCamera->nearDistance.setValue(1);
  fSoCamera->farDistance.setValue(100);
  fSoCamera->focalDistance.setValue(10);
  fSoSelection->addChild(fSoCamera);

  //fSoCameraSensor = new SoNodeSensor(CameraSensorCB,this);
  //fSoCameraSensor->attach(fSoCamera);

  fSoSelection->addChild(fG4OpenInventorSceneHandler.fRoot);

  // SoImageWriter should be the last.
  fSoImageWriter = new SoImageWriter();
  fSoImageWriter->fileName.setValue("g4out.ps");
  fSoSelection->addChild(fSoImageWriter);
}

G4OpenInventorViewer::~G4OpenInventorViewer () {
  //fSoCameraSensor->detach();
  //delete fSoCameraSensor;
  fSoSelection->unref();
}

void G4OpenInventorViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  if (
      //??fG4OpenInventorSceneHandler.fPODLList.size() == 0 ||
      // We need a test for empty scene graph, such as
      // staticRoot.size() or something??????????  See temporary fix
      // in contructor.  (John Allison Aug 2001)
      CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();
  }      
  fLastVP = fVP;
}
 
G4bool G4OpenInventorViewer::CompareForKernelVisit(G4ViewParameters& vp) {

  if (
      (vp.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (vp.GetRepStyle ()        != fVP.GetRepStyle ())        ||
      (vp.IsCulling ()          != fVP.IsCulling ())          ||
      (vp.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (vp.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (vp.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (vp.IsSection ()          != fVP.IsSection ())          ||
      // No need to visit kernel if section plane changes.
      (vp.IsCutaway ()          != fVP.IsCutaway ())          ||
      (vp.GetCutawayPlanes ().size () !=
                                 fVP.GetCutawayPlanes ().size ()) ||
      // No need to visit kernel if cutaway planes change.
      (vp.IsExplode ()          != fVP.IsExplode ())          ||
      (vp.GetNoOfSides ()       != fVP.GetNoOfSides ())
      ) {
      return true;;
  }
  if (vp.IsDensityCulling () &&
      (vp.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  if (vp.IsExplode () &&
      (vp.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;
      
  return false;
}

void G4OpenInventorViewer::ClearView () {
}

void G4OpenInventorViewer::SetView () {

  // Get G4 camera infos :
  const G4Point3D target
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const G4Vector3D& direction = fVP.GetViewpointDirection().unit();
  const G4Point3D cameraPosition = target + cameraDistance * direction;
  const G4double pnear = fVP.GetNearDistance (cameraDistance, radius);
  const G4double pfar  = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  //const G4Normal3D& up = fVP.GetUpVector ();  

  /*FIXME : G4OpenGLView::SetView does the below. Do we need that with IV ?
  G4Point3D gltarget;
  if (cameraDistance > 1.e-6 * radius) {
    gltarget = target;
  } else {
    gltarget = target - radius * fVP.GetViewpointDirection().unit();
  }
  */  

/*
  printf("debug : target : %g %g %g\n",target.x(),
                                       target.y(),
                                       target.z());
  printf("debug : dir : %g %g %g\n",direction.x(),
                                    direction.y(),
                                    direction.z());
  printf("debug : pos : %g %g %g\n",cameraPosition.x(),
                                    cameraPosition.y(),
                                    cameraPosition.z());
  printf("debug : near %g far %g\n",pnear,pfar);
*/

  // fSoCamera setup :
  fSoCamera->position.setValue((float)cameraPosition.x(),
                               (float)cameraPosition.y(),
                               (float)cameraPosition.z());

  //SbVec3f upFrom(0,1,0);
  //SbVec3f upTo((float)up.x(),(float)up.y(),(float)up.z());
  //fSoCamera->orientation.setValue(SbRotation(upFrom,upTo));
  SbVec3f sbTarget((float)target.x(),
                   (float)target.y(),
                   (float)target.z());
  sbTarget.normalize();
  // PointAt by keeping a up along y.
  fSoCamera->pointAt(sbTarget);

  //fSoCamera->height.setValue(10);
  fSoCamera->nearDistance.setValue((float)pnear);
  fSoCamera->farDistance.setValue((float)pfar);
  //fSoCamera->focalDistance.setValue((float)cameraDistance);
}

void G4OpenInventorViewer::DrawView () {
  //G4cout << "debug Iv::DrawViewer " <<G4endl;
  KernelVisitDecision();
  ProcessView();
  FinishView();
}

void G4OpenInventorViewer::ShowView () {
  fInteractorManager -> SecondaryLoop ();
}

/*
void G4OpenInventorViewer::CameraSensorCB(void* aThis,SoSensor*) { 
  G4OpenInventorViewer* This = (G4OpenInventorViewer*)aThis;

  SbVec3f direction;
  This->fSoCamera->orientation.getValue().multVec(SbVec3f(0,0,-1),direction);
  This->fVP.SetViewpointDirection
    (G4Vector3D(-direction[0],-direction[1],-direction[2]));

  SbVec3f pos = This->fSoCamera->position.getValue();
  SbVec3f target = pos + direction * This->fSoCamera->focalDistance.getValue();

  This->fVP.SetCurrentTargetPoint(G4Point3D(target[0],target[1],target[2]));
}
*/

void G4OpenInventorViewer::SelectionCB(
 void* aThis
,SoPath* aPath
) 
{
  G4OpenInventorViewer* This = (G4OpenInventorViewer*)aThis;
  SoNode* node = ((SoFullPath*)aPath)->getTail();
  G4String name((char*)node->getName().getString());
  G4String cls((char*)node->getTypeId().getName().getString());
  G4cout << "SoNode : " << node 
         << " SoType : " << cls 
         << " name : " << name 
         << G4endl;
/*FIXME : to explore
  if(node->isOfType(Geant4_SoPolyhedron::getClassTypeId())) {
    Geant4_SoPolyhedron* polyhedron = (Geant4_SoPolyhedron*)node;
    if(polyhedron->solid.getValue()==FALSE)
      polyhedron->solid.setValue(TRUE);
    else
      polyhedron->solid.setValue(FALSE);
  }*/
  This->fSoSelection->deselectAll();
}
/*
void G4OpenInventorViewer::DeselectionCB(
 void* aThis
,SoPath* aPath
) 
{
  //G4OpenInventorViewer* This = (G4OpenInventorViewer*)aThis;
  G4String name((char*)aPath->getTail()->getTypeId().getName().getString());
  G4cout << "Deselect : " << name << G4endl;
}
*/
//////////////////////////////////////////////////////////////////////////////
/// Menu items callbacks /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void G4OpenInventorViewer::Escape(){
  G4cout << "Escape..." <<G4endl;
  fInteractorManager->RequireExitSecondaryLoop (OIV_EXIT_CODE);
}

void G4OpenInventorViewer::WritePostScript(const G4String& aFile) {
  if(!fGL2PSAction) return;
  fGL2PSAction->setFileName(aFile.c_str());
  G4cout << "Produce " << aFile << "..." << G4endl;
  fGL2PSAction->enableFileWriting();
  ViewerRender();
  fGL2PSAction->disableFileWriting();
}

void G4OpenInventorViewer::WritePixmapPostScript(const G4String& aFile) {
  fSoImageWriter->fileName.setValue(aFile.c_str());
  //imageWriter->format.setValue(SoImageWriter::POST_SCRIPT);
  fSoImageWriter->enable();
  ViewerRender();
  fSoImageWriter->disable();
  if(fSoImageWriter->getStatus()) {
    G4cout << G4String(fSoImageWriter->fileName.getValue().getString()) 
           << " produced."
           << G4endl;
  } else {
    G4cout << G4String(fSoImageWriter->fileName.getValue().getString()) 
           << " not produced."
           << G4endl;
  }
}  

void G4OpenInventorViewer::WriteInventor(const G4String& aFile) {
  G4cout << "Produce " << aFile << "..." << G4endl;
  SoWriteAction writeAction;
  writeAction.getOutput()->openFile(aFile.c_str());
  writeAction.getOutput()->setBinary(FALSE);
  writeAction.apply(fSoSelection);
  writeAction.getOutput()->closeFile();
}

struct Counter {
 int fTriangles;
 int fLineSegments;
 int fPoints;
};

static void CountTrianglesCB(
 void* userData
,SoCallbackAction*
,const SoPrimitiveVertex*
,const SoPrimitiveVertex*,
const SoPrimitiveVertex*)
{
  Counter* counter = (Counter*)userData;
  counter->fTriangles++;
}

static void CountLineSegmentsCB(
 void* userData
,SoCallbackAction*
,const SoPrimitiveVertex*
,const SoPrimitiveVertex*)
{
  Counter* counter = (Counter*)userData;
  counter->fLineSegments++;
}

static void CountPointsCB(
 void* userData
,SoCallbackAction*
,const SoPrimitiveVertex*)
{
  Counter* counter = (Counter*)userData;
  counter->fPoints++;
}

void G4OpenInventorViewer::SceneGraphStatistics() {
  Counter counter;
  counter.fTriangles = 0;
  counter.fLineSegments = 0;
  counter.fPoints = 0;

  SoCallbackAction callbackAction;
  callbackAction.addTriangleCallback
    (SoShape::getClassTypeId(),CountTrianglesCB,(void*)&counter);
  callbackAction.addLineSegmentCallback
    (SoShape::getClassTypeId(),CountLineSegmentsCB,(void*)&counter);
  callbackAction.addPointCallback
    (SoShape::getClassTypeId(),CountPointsCB,(void*)&counter);
  callbackAction.apply(fSoSelection);

  SoCounterAction counterAction;
  counterAction.apply(fSoSelection);
  int nodes = counterAction.getCount();

  counterAction.setLookFor(SoCounterAction::TYPE);
  counterAction.setType(SoShape::getClassTypeId());
  counterAction.apply(fSoSelection);
  int shapes = counterAction.getCount();

  G4cout << "Number of triangles : " << counter.fTriangles << G4endl;
  G4cout << "Number of line segments : " << counter.fLineSegments << G4endl;
  G4cout << "Number of points : " << counter.fPoints << G4endl;
  G4cout << "Number of nodes : " << nodes << G4endl;
  G4cout << "Number of shapes : " << shapes << G4endl;
}

void G4OpenInventorViewer::EraseDetector() {
  fG4OpenInventorSceneHandler.fDetectorRoot->removeAllChildren();
}
void G4OpenInventorViewer::EraseEvent() {
  fG4OpenInventorSceneHandler.fTransientRoot->removeAllChildren();
}

void G4OpenInventorViewer::SetPreviewAndFull() {
  fG4OpenInventorSceneHandler.fPreviewAndFull = true;

  NeedKernelVisit();
  // DrawView does a ClearStore. Do not clear the transient store :
  SoSeparator* tmp = fG4OpenInventorSceneHandler.fTransientRoot;
  fG4OpenInventorSceneHandler.fTransientRoot = new SoSeparator;
  DrawView();  
  fG4OpenInventorSceneHandler.fTransientRoot->unref();
  fG4OpenInventorSceneHandler.fTransientRoot = tmp;
}

void G4OpenInventorViewer::SetPreview() {
  fG4OpenInventorSceneHandler.fPreviewAndFull = false;

  NeedKernelVisit();
  // DrawView does a ClearStore. Do not clear the transient store :
  SoSeparator* tmp = fG4OpenInventorSceneHandler.fTransientRoot;
  fG4OpenInventorSceneHandler.fTransientRoot = new SoSeparator;
  DrawView();  
  fG4OpenInventorSceneHandler.fTransientRoot->unref();
  fG4OpenInventorSceneHandler.fTransientRoot = tmp;
}

void G4OpenInventorViewer::SetSolid() {
  G4ViewParameters vp = GetViewParameters();
  G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
  //From G4VisCommandsViewerSet : /vis/viewer/set/style solid.
  switch (existingStyle) {
  case G4ViewParameters::wireframe:
    vp.SetDrawingStyle(G4ViewParameters::hsr);
    break;
  case G4ViewParameters::hlr:
    vp.SetDrawingStyle(G4ViewParameters::hlhsr);
    break;
  case G4ViewParameters::hsr:
    break;
  case G4ViewParameters::hlhsr:
    break;
  }
  SetViewParameters(vp);
/*FIXME : fSoCameraSensor : when ready, should use the below.
  if (vp.IsAutoRefresh()) {
    //G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    SetView();
    ClearView();
    //FIXME : do not clear transient store.
    DrawView();  
  }
*/
  //FIXME : waiting fSoCameraSensor.
  fModified = false;
  // DrawView does a ClearStore. Do not clear the transient store :
  SoSeparator* tmp = fG4OpenInventorSceneHandler.fTransientRoot;
  fG4OpenInventorSceneHandler.fTransientRoot = new SoSeparator;
  DrawView();  
  fG4OpenInventorSceneHandler.fTransientRoot->unref();
  fG4OpenInventorSceneHandler.fTransientRoot = tmp;
}
void G4OpenInventorViewer::SetWireFrame() {
  G4ViewParameters vp = GetViewParameters();
  G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
  //From G4VisCommandsViewerSet : /vis/viewer/set/style wire.
  switch (existingStyle) {
  case G4ViewParameters::wireframe:
    break;
  case G4ViewParameters::hlr:
    break;
  case G4ViewParameters::hsr:
    vp.SetDrawingStyle(G4ViewParameters::wireframe);
    break;
  case G4ViewParameters::hlhsr:
    vp.SetDrawingStyle(G4ViewParameters::hlr);
    break;
  }
  SetViewParameters(vp);
/*FIXME : fSoCameraSensor : when ready, should use the below.
  if (vp.IsAutoRefresh()) {
    //G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    SetView();
    ClearView();
    //FIXME : do not clear transient store.
    DrawView();
  }
*/
  //FIXME : waiting fSoCameraSensor.
  fModified = false;
  // DrawView does a ClearStore. Do not clear the transient store :
  SoSeparator* tmp = fG4OpenInventorSceneHandler.fTransientRoot;
  fG4OpenInventorSceneHandler.fTransientRoot = new SoSeparator;
  DrawView();  
  fG4OpenInventorSceneHandler.fTransientRoot->unref();
  fG4OpenInventorSceneHandler.fTransientRoot = tmp;
}
void G4OpenInventorViewer::UpdateScene() {
  fG4OpenInventorSceneHandler.ClearStore();
  //FIXME : fSoCameraSensor : when ready, should use the below.
  //SetView();
  ClearView();
  DrawView();
  ShowView();
}
G4String G4OpenInventorViewer::Help(const G4String& aTopic) {
  if(aTopic=="controls") {
    return G4String("\
Controls on an Inventor examiner viewer are :\n\
- in picking mode (cursor is the upper left arrow)\n\
  Ctrl + pick a volume : see daughters.\n\
  Shift + pick a volume : see mother.\n\
- in viewing mode (cursor is the hand)\n\
  Left-button + pointer move : rotate.\n\
  Ctrl+Left-button + pointer move : pan.\n\
  Ctrl+Shift+Left-button + pointer move : scale.\n\
  Middle-button + pointer move : pan.\n\
  Right-button : popup menu.\n");
  } else {
    return "";
  }
}

#endif



