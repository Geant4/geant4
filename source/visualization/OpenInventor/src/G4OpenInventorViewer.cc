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
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoWriteAction.h>

#include "HEPVis/nodes/SoImageWriter.h"
#include "HEPVis/actions/SoGL2PSAction.h"

#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4VInteractorManager.hh"

G4OpenInventorViewer::G4OpenInventorViewer(
 G4OpenInventorSceneHandler& sceneHandler
,const G4String& name)
:G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
,fG4OpenInventorSceneHandler(sceneHandler)
,fInteractorManager(0)
,fSoSelection(0)
,fSoImageWriter(0)
,fGL2PSAction(0) //To be set be suclass.
{
  fNeedKernelVisit = true;  //?? Temporary, until KernelVisitDecision fixed.

  // For viewing of all objects by default :
  fDefaultVP.SetCulling(false);
  fVP.SetCulling(false);

  fInteractorManager = 
    ((G4OpenInventor*)fG4OpenInventorSceneHandler.GetGraphicsSystem())->
    GetInteractorManager();

  // Main user scene graph root sent to the viewers.
  fSoSelection = new SoSelection;
  fSoSelection->policy = SoSelection::SINGLE;
  fSoSelection->ref();
  fSoSelection->addChild(fG4OpenInventorSceneHandler.fRoot);

  fSoImageWriter = new SoImageWriter();
  fSoImageWriter->fileName.setValue("g4out.ps");
  fSoSelection->addChild(fSoImageWriter);
}

G4OpenInventorViewer::~G4OpenInventorViewer () {
  if(fSoSelection) fSoSelection->unref();
}

void G4OpenInventorViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  if (
      //??fG4OpenInventorSceneHandler.fPODLList.size() == 0 ||
      // We need a test for empty scene graph, such as
      // staticRoot.size() or something??????????  See temporary fix
      // in contructor.  (John Allison Aug 2001)
      CompareForKernelVisit(fG4OpenInventorSceneHandler.fLastVP)  ||
      CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();
  }      
  fLastVP = fVP;
  fG4OpenInventorSceneHandler.fLastVP = fVP;
}
 
G4bool G4OpenInventorViewer::CompareForKernelVisit(G4ViewParameters&) {

  if (
      (fLastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (fLastVP.GetRepStyle ()        != fVP.GetRepStyle ())        ||
      (fLastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (fLastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (fLastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (fLastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (fLastVP.IsSection ()          != fVP.IsSection ())          ||
      // No need to visit kernel if section plane changes.
      (fLastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      (fLastVP.GetCutawayPlanes ().size () !=
                                 fVP.GetCutawayPlanes ().size ()) ||
      // No need to visit kernel if cutaway planes change.
      (fLastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (fLastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())
      ) {
      return true;;
  }
  if (fLastVP.IsDensityCulling () &&
      (fLastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  if (fLastVP.IsExplode () &&
      (fLastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;
      
  return false;
}

void G4OpenInventorViewer::ClearView () {
}

void G4OpenInventorViewer::SetView () {
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

static void CountTriangleCB(
 void* userData
,SoCallbackAction*
,const SoPrimitiveVertex*
,const SoPrimitiveVertex*
,const SoPrimitiveVertex*)
{
  int* number = (int*)userData;
  (*number)++;
}

void G4OpenInventorViewer::CountTriangles() {
  SoCallbackAction callbackAction;
  int trianglen = 0;
  fSoSelection->ref();
  callbackAction.addTriangleCallback
    (SoShape::getClassTypeId(),CountTriangleCB,(void*)&trianglen);
  callbackAction.apply(fSoSelection);
  fSoSelection->unrefNoDelete();
  G4cout << "Number of triangles : " << trianglen << G4endl;
}

void G4OpenInventorViewer::EraseDetector() {
  fG4OpenInventorSceneHandler.fDetectorRoot->removeAllChildren();
}
void G4OpenInventorViewer::EraseEvent() {
  fG4OpenInventorSceneHandler.fTransientRoot->removeAllChildren();
}

#endif
