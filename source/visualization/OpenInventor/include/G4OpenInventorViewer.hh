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

#ifndef G4OPENINVENTORVIEWER_HH
#define G4OPENINVENTORVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4VViewer.hh"

class SoSelection;
class SoPath;
class SoCamera;
class SoSensor;
class SoNodeSensor;

class Geant4_SoImageWriter;
class Geant4_SoGL2PSAction;
class G4OpenInventorSceneHandler;
class G4VInteractorManager;
class SbVec3f;

//
// Base class for various OpenInventorView classes.
//
class G4OpenInventorViewer: public G4VViewer {
public: //G4VViewer
  virtual void DrawView();
  virtual void ShowView();
  virtual void ClearView();
  virtual void SetView();
  virtual void KernelVisitDecision();
public:
  G4OpenInventorViewer(G4OpenInventorSceneHandler& scene,
		       const G4String& name = "");
  virtual ~G4OpenInventorViewer();
protected:
  virtual void ViewerRender() = 0;
  virtual SoCamera* GetCamera() = 0;
  void Escape();
  void WritePostScript(const G4String& file = "g4out.ps");
  void WritePDF(const G4String& file = "g4out.pdf");
  void WritePixmapPostScript(const G4String& file = "g4out.ps");
  void WriteInventor(const G4String& file = "g4out.iv");
  void SceneGraphStatistics();
  void EraseDetector();
  void EraseEvent();
  void SetPreviewAndFull();
  void SetPreview();
  void SetSolid();
  void SetWireFrame();
  void SetReducedWireFrame(bool);
  void UpdateScene();
  G4String Help(const G4String& topic = "controls");
private:
  //static void SelectionCB(void*,SoPath*);
  //static void DeselectionCB(void*,SoPath*);
  static void GroupCameraSensorCB(void*,SoSensor*);
  static void CameraSensorCB(void*,SoSensor*);
  static void pointAt(SoCamera*,const SbVec3f & targetpoint, const SbVec3f & upvector);
  static void lookAt(SoCamera*,const SbVec3f & dir, const SbVec3f & up);
  static void lookedAt(SoCamera*,SbVec3f & dir, SbVec3f & up);
private:
  G4bool CompareForKernelVisit(G4ViewParameters&);
  void DrawDetector();
private:
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
protected:
  static void SelectionCB(void*,SoPath*);
  G4OpenInventorSceneHandler& fG4OpenInventorSceneHandler;
  G4VInteractorManager* fInteractorManager;
  SoSelection* fSoSelection;
  Geant4_SoImageWriter* fSoImageWriter;
  Geant4_SoGL2PSAction* fGL2PSAction;
  SoNodeSensor* fGroupCameraSensor;
  SoNodeSensor* fCameraSensor;
};

#endif

#endif
