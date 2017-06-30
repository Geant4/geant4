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
// $Id: G4OpenInventorViewer.cc 102233 2017-01-13 15:59:06Z gcosmo $

#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "G4OpenInventorViewer.hh"

#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoShape.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/sensors/SoNodeSensor.h>

#include "HEPVis/nodes/SoImageWriter.h"
#include "HEPVis/actions/SoGL2PSAction.h"
#include "HEPVis/actions/SoCounterAction.h"
#include "HEPVis/actions/SoAlternateRepAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"
#include "G4Scene.hh"
#include "Geant4_SoPolyhedron.h"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"
#include "G4AttHolder.hh"

G4OpenInventorViewer::G4OpenInventorViewer(
 G4OpenInventorSceneHandler& sceneHandler
,const G4String& name)
:G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
,fG4OpenInventorSceneHandler(sceneHandler)
,fInteractorManager(0)
,fSoSelection(0)
,fSoImageWriter(0)
,fGL2PSAction(0) //To be set be suclass.
,fGroupCameraSensor(0)
,fCameraSensor(0)
{
  fNeedKernelVisit = true;  //?? Temporary, until KernelVisitDecision fixed.

  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);
  fVP.SetPicking(true);
  fDefaultVP.SetPicking(true);

  //FIXME : G.Barrand : not convinced that we have to rm culling.
  // For viewing of all objects by default :
  //fDefaultVP.SetCulling(false);
  //fVP.SetCulling(false);

  fInteractorManager = 
    ((G4OpenInventor*)fG4OpenInventorSceneHandler.GetGraphicsSystem())->
    GetInteractorManager();

  // Main user scene graph root sent to the viewers.
  fSoSelection = new SoSelection;
  fSoSelection->ref();
  fSoSelection->addSelectionCallback(SelectionCB,this);
  //fSoSelection->addDeselectionCallback(DeselectionCB,this);
  fSoSelection->policy = SoSelection::SINGLE;

  SoGroup* group = new SoGroup;
  fSoSelection->addChild(group);

  //  Have a camera under fSoSelection in order
  // that the below SceneGraphSensor be notifed
  // when the viewer changes the camera type.
  //  But we put the camera under a SoGroup so that
  // the SceneGraphSensor be not triggered at each change
  // under the fG4OpenInventorSceneHandler.fRoot.
  SoOrthographicCamera* camera = new SoOrthographicCamera;
  camera->viewportMapping.setValue(SoCamera::ADJUST_CAMERA);
  //camera->aspectRatio.setValue(10);
  camera->position.setValue(0,0,10);
  camera->orientation.setValue(SbRotation(SbVec3f(0,1,0),0));
  camera->height.setValue(10);
  camera->nearDistance.setValue(1);
  camera->farDistance.setValue(100);
  camera->focalDistance.setValue(10);
  group->addChild(camera);

 {SoInput soInput;
   if(soInput.openFile("g4view.iv",TRUE)) {
    SoSeparator* separator = SoDB::readAll(&soInput);
    if(separator) fSoSelection->addChild(separator);
  }}

  fSoSelection->addChild(fG4OpenInventorSceneHandler.fRoot);

  // SoImageWriter should be the last.
  fSoImageWriter = new SoImageWriter();
  fSoImageWriter->fileName.setValue("g4out.ps");
  fSoSelection->addChild(fSoImageWriter);

  // Sensors :
  // To detect that the viewer had changed the camera type :
  fGroupCameraSensor = new SoNodeSensor(GroupCameraSensorCB,this);
  fGroupCameraSensor->setPriority(0);//Needed in order to do getTriggerNode()
  fGroupCameraSensor->attach(group);

  fCameraSensor = new SoNodeSensor(CameraSensorCB,this);
  fCameraSensor->setPriority(0);//Needed in order to do getTriggerNode()
}

G4OpenInventorViewer::~G4OpenInventorViewer () {
  fCameraSensor->detach();
  delete fCameraSensor;
  fGroupCameraSensor->detach();
  delete fGroupCameraSensor;
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
      (vp.IsAuxEdgeVisible ()   != fVP.IsAuxEdgeVisible ())   ||
      (vp.IsCulling ()          != fVP.IsCulling ())          ||
      (vp.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (vp.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (vp.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (vp.IsSection ()          != fVP.IsSection ())          ||
      (vp.IsCutaway ()          != fVP.IsCutaway ())          ||
      // This assumes use of generic clipping (sectioning, slicing,
      // DCUT, cutaway).  If a decision is made to implement locally,
      // this will need changing.  See G4OpenGLViewer::SetView,
      // G4OpenGLStoredViewer.cc::CompareForKernelVisit and
      // G4OpenGLStoredSceneHander::CreateSection/CutawayPolyhedron.
      (vp.IsExplode ()          != fVP.IsExplode ())          ||
      (vp.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
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
      // If G4OpenInventor ever introduces VAMs, the following might need
      // changing to a complete comparison, i.e., remove ".size()".  See
      // G4OpenGLStoredViewer::CompareForKernelVisit.
      (vp.GetVisAttributesModifiers().size() !=
       fVP.GetVisAttributesModifiers().size())
      )
    return true;

  if (vp.IsDensityCulling () &&
      (vp.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  if (vp.IsSection () &&
      (vp.GetSectionPlane () != fVP.GetSectionPlane ()))
    return true;

  if (vp.IsCutaway ()) {
    if (vp.GetCutawayPlanes ().size () !=
	fVP.GetCutawayPlanes ().size ()) return true;
    for (size_t i = 0; i < vp.GetCutawayPlanes().size(); ++i)
      if (vp.GetCutawayPlanes()[i] != fVP.GetCutawayPlanes()[i])
	return true;
  }

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
  //const G4double pnear = fVP.GetNearDistance (cameraDistance, radius);
  //const G4double pfar  = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  const G4Normal3D& up = fVP.GetUpVector ();  

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
  //printf("debug : near %g far %g\n",pnear,pfar);
*/

  SoCamera* camera = GetCamera();
  if(!camera) return;

  // viewer camera setup :
  camera->position.setValue((float)cameraPosition.x(),
                               (float)cameraPosition.y(),
                               (float)cameraPosition.z());

  SbVec3f sbTarget((float)target.x(),
                   (float)target.y(),
                   (float)target.z());
  SbVec3f sbUp((float)up.x(),
	       (float)up.y(),
	       (float)up.z());
  sbUp.normalize();
  // Need Coin's camera->pointAt(sbTarget,sbUp); not in the SGI API
  // Stole Coin's code...
  pointAt(camera,sbTarget,sbUp);

  //camera->height.setValue(10);
  //camera->nearDistance.setValue((float)pnear);
  //camera->farDistance.setValue((float)pfar);
  //camera->focalDistance.setValue((float)cameraDistance);

  if(camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
    if (fVP.GetFieldHalfAngle() == 0.) {
      //FIXME : ((SoOrthographicCamera*)camera)->height.setValue();
      //FIXME : (Don't think we have to do that.)
    } else {
      //FIXME : Have to set a perspective camera !
      //FIXME : viewer->setCameraType(SoPerspectiveCamera::getClassTypeId())
      //FIXME : ((SoPerspectiveCamera*)camera)->heightAngle.setValue
      //FIXME :   (2.*fVP.GetFieldHalfAngle());
    }
  } else if(camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
    if (fVP.GetFieldHalfAngle() == 0.) {
      //FIXME : Have to set an orthographic camera !
      //FIXME : viewer->setCameraType(SoOrthographicCamera::getClassTypeId())
    } else {
      //FIXME : ((SoPerspectiveCamera*)camera)->heightAngle.setValue
      //FIXME :   (2.*fVP.GetFieldHalfAngle());
    }
  }
}

//COIN_FUNCTION_EXTENSION
void
G4OpenInventorViewer::pointAt(SoCamera* camera,const SbVec3f & targetpoint, const SbVec3f & upvector)
{
  SbVec3f dir = targetpoint - camera->position.getValue();
  if (dir.normalize() == 0.0f) return;
  lookAt(camera,dir, upvector);
}

//COIN_FUNCTION
// Private method that calculates a new orientation based on camera
// direction and camera up vector. Vectors must be unit length.
void
G4OpenInventorViewer::lookAt(SoCamera* camera,const SbVec3f & dir, const SbVec3f & up)
{
  SbVec3f z = -dir;
  SbVec3f y = up;
  SbVec3f x = y.cross(z);

  // recompute y to create a valid coordinate system
  y = z.cross(x);

  // normalize x and y to create an orthonormal coord system
  y.normalize();
  x.normalize();

  // create a rotation matrix
  SbMatrix rot = SbMatrix::identity();
  rot[0][0] = x[0];
  rot[0][1] = x[1];
  rot[0][2] = x[2];

  rot[1][0] = y[0];
  rot[1][1] = y[1];
  rot[1][2] = y[2];

  rot[2][0] = z[0];
  rot[2][1] = z[1];
  rot[2][2] = z[2];

  camera->orientation.setValue(SbRotation(rot));
}

void
G4OpenInventorViewer::lookedAt(SoCamera* camera,SbVec3f & dir, SbVec3f & up)
{
  SbRotation rot = camera->orientation.getValue();
  SbMatrix mrot; rot.getValue(mrot);

  SbVec3f x, y, z;

  // create a rotation matrix
  x[0] = mrot[0][0];
  x[1] = mrot[0][1];
  x[2] = mrot[0][2];

  y[0] = mrot[1][0];
  y[1] = mrot[1][1];
  y[2] = mrot[1][2];

  z[0] = mrot[2][0];
  z[1] = mrot[2][1];
  z[2] = mrot[2][2];

  dir = -z;
  dir.normalize();
  up = SbVec3f(0.f,1.f,0.f);  // Choose y-axis if possible.
  if (std::abs(up.dot(z)) > 1.e-6) {
    up = y;
    up.normalize();
  }
}

void G4OpenInventorViewer::DrawView () {
  //G4cout << "debug Iv::DrawViewer " <<G4endl;
  if (!fNeedKernelVisit) KernelVisitDecision();
  ProcessView();
  FinishView();
}

void G4OpenInventorViewer::ShowView () {
  fInteractorManager -> SecondaryLoop ();
}

void G4OpenInventorViewer::GroupCameraSensorCB(void* aThis,SoSensor* aSensor){ 
  G4OpenInventorViewer* This = (G4OpenInventorViewer*)aThis;

  SoNode* node = ((SoNodeSensor*)aSensor)->getTriggerNode();
  //printf("debug : GroupCameraSensorCB %s\n",
  //node->getTypeId().getName().getString());

  if(node->isOfType(SoCamera::getClassTypeId())) {
    // Viewer had changed the camera type, 
    // attach the fCameraSensor to the new camera.
    SoCamera* camera = (SoCamera*)node;
    This->fCameraSensor->detach();
    This->fCameraSensor->attach(camera);
  }

}

void G4OpenInventorViewer::CameraSensorCB(void* aThis,SoSensor* aSensor) { 
  G4OpenInventorViewer* This = (G4OpenInventorViewer*)aThis;

  //printf("debug : CameraSensorCB\n");

  SoNode* node = ((SoNodeSensor*)aSensor)->getTriggerNode();

  if(node->isOfType(SoCamera::getClassTypeId())) {
    SoCamera* camera = (SoCamera*)node;

    SbVec3f direction, up;
    lookedAt(camera,direction, up);
    This->fVP.SetViewpointDirection
      (G4Vector3D(-direction[0],-direction[1],-direction[2]));
    This->fVP.SetUpVector(G4Vector3D(up[0],up[1],up[2]));

    SbVec3f pos = camera->position.getValue();
    SbVec3f target = pos + direction * camera->focalDistance.getValue();

    This->fVP.SetCurrentTargetPoint(G4Point3D(target[0],target[1],target[2]));
  }
}

void G4OpenInventorViewer::SelectionCB(
 void* aThis
,SoPath* aPath
) 
{
  G4OpenInventorViewer* This = (G4OpenInventorViewer*)aThis;
  SoNode* node = ((SoFullPath*)aPath)->getTail();
  G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
  if(attHolder && attHolder->GetAttDefs().size()) {
    for (size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
      G4cout << G4AttCheck(attHolder->GetAttValues()[i],
			   attHolder->GetAttDefs()[i]);
    }
  } else {
    G4String name((char*)node->getName().getString());
    G4String cls((char*)node->getTypeId().getName().getString());
    G4cout << "SoNode : " << node 
	   << " SoType : " << cls 
	   << " name : " << name 
	   << G4endl;
    G4cout << "No attributes attached." << G4endl;
  }
  /*FIXME : to explore (need different button - this is used for picking.
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

void G4OpenInventorViewer::DrawDetector() {
  /* Replace this... - JA
  // DrawView does a ClearStore. Do not clear the transient store :
  SoSeparator* tmp = fG4OpenInventorSceneHandler.fTransientRoot;
  fG4OpenInventorSceneHandler.fTransientRoot = new SoSeparator;
  if (!fNeedKernelVisit) KernelVisitDecision();
  ProcessView();
  fG4OpenInventorSceneHandler.fTransientRoot->unref();
  fG4OpenInventorSceneHandler.fTransientRoot = tmp;
  */
  // ...by this... - JA
  DrawView();
}

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
  fGL2PSAction->setExportImageFormat(GL2PS_EPS);
  // Use gl2ps default buffer (2048*2048)
  fGL2PSAction->setBufferSize(0);
  G4cout << "Produce " << aFile << "..." << G4endl;
  if (fGL2PSAction->enableFileWriting()) {
    ViewerRender();
    fGL2PSAction->disableFileWriting();
  }
  fGL2PSAction->resetBufferSizeParameters();
}

void G4OpenInventorViewer::WritePDF(const G4String& aFile) {
  if(!fGL2PSAction) return;
  fGL2PSAction->setFileName(aFile.c_str());
  fGL2PSAction->setExportImageFormat(GL2PS_PDF);
  // Use gl2ps default buffer (2048*2048)
  fGL2PSAction->setBufferSize(0);
  G4cout << "Produce " << aFile << "..." << G4endl;
  if (fGL2PSAction->enableFileWriting()) {
    ViewerRender();
    fGL2PSAction->disableFileWriting();
  }
  fGL2PSAction->resetBufferSizeParameters();
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

  SbBool genAlternateRep = TRUE;
  //SbBool binary = FALSE;
  SbBool binary = TRUE;
  SoAlternateRepAction alternateRepAction;
  if(genAlternateRep==TRUE) {
    alternateRepAction.setGenerate(TRUE); //Clear alternate reps.
    alternateRepAction.apply(fSoSelection);
  }

  SoWriteAction writeAction;
  writeAction.getOutput()->openFile(aFile.c_str());
  writeAction.getOutput()->setBinary(binary);
  writeAction.apply(fSoSelection);
  writeAction.getOutput()->closeFile();

  if(genAlternateRep==TRUE) {
    alternateRepAction.setGenerate(FALSE); //Clear alternate reps.
    alternateRepAction.apply(fSoSelection);
  }



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
  DrawDetector();
}

void G4OpenInventorViewer::SetPreview() {
  fG4OpenInventorSceneHandler.fPreviewAndFull = false;

  NeedKernelVisit();
  DrawDetector();
}

// When ViewParameter <-> SoCamera mapping ready 
// uncomment the below
//#define USE_SET_VIEW

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
  DrawDetector();
}
void G4OpenInventorViewer::SetWireFrame() {
  G4ViewParameters vp = GetViewParameters();
  G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
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
  DrawDetector();
}


void G4OpenInventorViewer::SetReducedWireFrame(bool aValue) {
  G4ViewParameters vp = GetViewParameters();

  // Set the wire frame kind :
  vp.SetAuxEdgeVisible(!aValue);

  // Set wire frame :
  G4ViewParameters::DrawingStyle existingStyle = vp.GetDrawingStyle();
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
  NeedKernelVisit(); // Just in case it was alread in wire framw.
  DrawDetector();
}

void G4OpenInventorViewer::UpdateScene() {
  /* Replace this... - JA
  fG4OpenInventorSceneHandler.ClearStore();
  ClearView();
  if (!fNeedKernelVisit) KernelVisitDecision();
  ProcessView();
  ShowView();
  */
  // ...by this - JA
  NeedKernelVisit();
  DrawView();
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
