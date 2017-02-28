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
// Open Inventor Xt Extended Viewer - 30 Oct 2012
// Rastislav Ondrasek, Pierre-Luc Gagnon, Frederick Jones TRIUMF

#ifndef HookEventProcState_H
#define HookEventProcState_H 1
#include "G4VStateDependent.hh"

class G4OpenInventorXtExaminerViewer;

class HookEventProcState : public G4VStateDependent
{
private:
   G4OpenInventorXtExaminerViewer *viewer;
public:
   HookEventProcState(G4OpenInventorXtExaminerViewer*);
   ~HookEventProcState();

   virtual G4bool Notify(G4ApplicationState requiredState);
};
#endif /* HookEventProcState_H */


#ifndef G4OPENINVENTORXTEXAMINERVIEWER_HH
#define G4OPENINVENTORXTEXAMINERVIEWER_HH

#include <map>
#include <vector>
#include <fstream>
#include <Inventor/SbLinear.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Inventor/events/SoKeyboardEvent.h>

class SoCoordinate3;
class SoFont;
class SoText2;
class SoPointSet;

class G4OpenInventorXtExaminerViewer : public SoXtExaminerViewer {

  friend class G4OpenInventorXtExaminerViewerMessenger;
  // FWJ
  friend class G4OpenInventorXtExtendedViewer;

private:
  Widget prevViewPtButton, nextViewPtButton;
  Widget menuBar, fileMenu, openFileDialog, newFileDialog,
  loadRefCoordsDialog, saveRefCoordsDialog,
  loadSceneGraphDialog, saveSceneGraphDialog,
  viewPtSelection, listsDialog, myShellDialog, myViewPtList, myElementList;

  static G4OpenInventorXtExaminerViewer *viewer;
  void (*escapeCallback)(void *);
  void * examinerObject;
  SbBool lshiftdown, rshiftdown, lctrldown, rctrldown;
  
public:

  // Same constructor as the ExaminerViewer
  G4OpenInventorXtExaminerViewer(Widget parent = NULL,
	     const char *name = NULL,
	     SbBool embed = TRUE, 
	     SoXtFullViewer::BuildFlag flag = BUILD_ALL,
	     SoXtViewer::Type type = BROWSER);

  ~G4OpenInventorXtExaminerViewer();

  template <class T> void parseString(T &t, const std::string &s, bool &error);

  Widget addMenu(std::string name);
  void addButton(Widget menu, std::string name, XtCallbackProc);
  Widget getMenuBar() { return menuBar; }
  Widget getMenu() { return fileMenu; }
  void warningMsgDialog(std::string, String, XtCallbackProc);
  bool warningFlag;

  std::string saveScenegraphFileName;
  Widget saveScenegraphWidget;
  std::string saveRefCoordsFileName;
  Widget saveRefCoordsWidget;

  Widget createScale(Widget, char *, int, float);
  void addEscapeCallback(void (*cb)(void *), void *);
  bool abbrOutputFlag;
  bool pickRefPathFlag;
  bool viewingBeforePickRef;
   // FWJ
   //  SoNode * superimposition;
		     
protected:
  // Same constructor as the ExaminerViewer 
  G4OpenInventorXtExaminerViewer(Widget parent,
	     const char *name,
	     SbBool embed,
	     SoXtFullViewer::BuildFlag flag,
	     SoXtViewer::Type type,
	     SbBool build);

  // Overloaded for adding the MenuBar
  Widget buildWidget(Widget parent);
  // Overloaded so additional buttons can be added
  virtual void createViewerButtons (Widget parent, SbPList * buttonlist);
  // Overloaded for catching various keyboard events  
  virtual SbBool processSoEvent(const SoEvent * const event);
  void moveCamera(float dist = 0, bool lookdown = false);
  std::string curEltName;
  SbVec3f camUpVec;
  SbVec3f camDir;
  void rotateCamera();
  void updateViewParams(SoKeyboardEvent::Key);
  bool loadViewPts();
  virtual void afterRealizeHook();

private:
  // Each constructor calls this generic constructor
  void constructor(const SbBool build);

  // FWJ DISABLED
  //  static G4OpenInventorXtExaminerViewer *getObject();

  HookEventProcState *hookBeamOn;
  friend class HookEventProcState;
  bool newEvents;
  static void sceneChangeCB(void *, SoSensor *);
    
  void setViewPt();
  void writeViewPtIdx();
  void cleanUpAfterPrevFile();
  
  void popUpFileSelDialog(Widget&, std::string, std::string, XtCallbackProc);
  static void cancelFileSelDialogCB(Widget, XtPointer, XtPointer);
  static void openViewPtFileCB(Widget, XtPointer, XtPointer);
  static void viewPtFileSelectedCB(Widget, XtPointer, XtPointer);
  static void newViewPtFileCB(Widget, XtPointer, XtPointer);
  static void createNewVPFileCB(Widget, XtPointer, XtPointer);
  static void overwriteFileCB(Widget, XtPointer, XtPointer);
  static void loadRefCoordsDialogCB(Widget, XtPointer, XtPointer);	//pop file dialog
  static void loadRefCoordsCB(Widget, XtPointer, XtPointer);			//execute loading
  static void saveRefCoordsDialogCB(Widget, XtPointer, XtPointer);	//pop file dialog
  static void saveRefCoordsCB(Widget, XtPointer, XtPointer);			//execute saving
  static void saveRefCoordsOverWriteCB(Widget, XtPointer, XtPointer);
  static void loadSceneGraphDialogCB(Widget, XtPointer, XtPointer);
  static void loadSceneGraphCB(Widget, XtPointer, XtPointer);
  static void saveSceneGraphDialogCB(Widget, XtPointer, XtPointer);
  static void saveSceneGraphCB(Widget, XtPointer, XtPointer);
  static void saveSceneGraphOverWriteCB(Widget, XtPointer, XtPointer);
  static void mouseoverCB(void *aThis, SoEventCallback *eventCB);
  static void pickingCB(void *aThis, SoEventCallback *eventCB);

  		     
  // Viewpoint operations
  void addViewPoints();
  static void closeListsDialogCB(Widget, XtPointer, XtPointer);
  static void loadBookmarkCB(Widget, XtPointer, XtPointer);
  static void renameBookmarkCB(Widget, XtPointer, XtPointer);
  void renameViewPt(char *vpName);
  static void sortBookmarksCB(Widget, XtPointer, XtPointer);
  void sortViewPts(std::vector<std::string>);
  static void deleteBookmarkCB(Widget, XtPointer, XtPointer);
  static void deleteViewPtCB(Widget, XtPointer, XtPointer);
  void deleteViewPt(char *vpName = NULL);
    
  // Animation
  static void animateRefParticleCB(Widget, XtPointer, XtPointer);
  static void animateSensorCB(void *, SoSensor *);
  static void animateSensorRotationCB(void *, SoSensor *);
  void animateRefParticle();
  void saveCurCamera();
  void restoreCamera();
  double animateBtwPtsPeriod, speedStep;
  void incSpeed();
  void decSpeed();
  
  SoTimerSensor *animateSensor;
  SoTimerSensor *animateSensorRotation;
  SoNodeSensor *sceneChangeSensor;
  SbVec3f camStartPos, camEndPos;
  SbRotation camStartOrient, camEndOrient;
  
  static void prevViewPtCB(Widget, XtPointer, XtPointer);
  static void nextViewPtCB(Widget, XtPointer, XtPointer);
  static void saveViewPtCB(Widget, XtPointer, XtPointer);
  static void abbrOutputCB(Widget, XtPointer, XtPointer);
  static void pickRefPathCB(Widget, XtPointer, XtPointer);
  static void switchWireFrameCB(Widget, XtPointer, XtPointer);
  static void constructListsDialog(Widget, XtPointer, XtPointer);
  void saveViewPt(char *name);
  

  static void lookAtSceneElementCB(Widget, XtPointer, XtPointer);
  static void cancelSceneElementSelectionCB(Widget, XtPointer, XtPointer);

  void setReferencePath(SoLineSet*, SoCoordinate3*, bool append = false);
  void setReferencePathZPos();
  void findAndSetRefPath();
  SoCoordinate3* getCoordsNode(SoFullPath *path);
  void getSceneElements(); // reads elements from the scene graph
  float sqrlen(const SbVec3f&);
  void distanceToTrajectory(const SbVec3f&, float&, SbVec3f&, int&);
  void sortElements();
  void createElementsList(Widget);
  static void closeMainWindowCB(Widget, XtPointer, XtPointer);
  void evenOutRefParticlePts();

  static void gotoRefPathStartCB(Widget, XtPointer, XtPointer);
  void gotoRefPathStart();
  static void invertRefPathCB(Widget, XtPointer, XtPointer);
  void invertRefPath();

  enum CameraType {
    PERSPECTIVE,
    ORTHOGRAPHIC
  };

    
  enum State {
    GENERAL,
    BEAMLINE,
    VIEWPOINT,
    ANIMATION,
    REVERSED_ANIMATION,
    PAUSED_ANIMATION,
    ROTATING
  };

  // For storing the view point
  struct viewPtData {
	  char* viewPtName;
	  int viewportMapping;
	  SbVec3f position;
	  SbRotation orientation;
	  float	aspectRatio;
	  float nearDistance;
	  float	farDistance;
	  float	focalDistance;
	  CameraType camType;
	  float	height;
  };

   // FWJ removed unneeded assignment operator
   struct sceneElement {
      std::string name;
      SoFullPath* path;
      SbVec3f center;
      float closestPointZCoord;
   };

   struct elementForSorting {
      float closestPointZCoord;
      SbVec3f closestPoint;
      float smallestDistance;
      float distanceToBeamlineStart;
      std::string name;

      bool operator<(elementForSorting const &other) const
      {
         if (closestPointZCoord < other.closestPointZCoord)
            return true;
         if (closestPointZCoord > other.closestPointZCoord)
            return false;

         // otherwise closestPointZCoord == other.closestPointZCoord.
         // Compare the distances from the center of the element to
         // the start of the beamline.
         if (distanceToBeamlineStart < other.distanceToBeamlineStart)
            return true;
         if (distanceToBeamlineStart > other.distanceToBeamlineStart)
            return false;

         // In case both closestPointZCoord and smallestDistance are
         // equal, we have two exactly overlapping elements, if so
         // the order doesn't matter
         return true;
      }
  };

  bool zcoordSetFlag;

  std::vector<sceneElement> sceneElements;
  std::vector<viewPtData> viewPtList;
  std::string fileName;
  std::ifstream fileIn;
  std::ofstream fileOut; 
  int viewPtIdx;
  int MAX_VP_IDX;
  int MAX_VP_NAME;
  
  // For storing coordinate points of the reference particle
  std::vector<SbVec3f> refParticleTrajectory;
  // For displaying distance during anim and beamline modes
  std::vector<float> refZPositions;

  int refParticleIdx;
  int prevRefIdx;
  float distance;
  State currentState, prevState, beforePausing;
  char* curViewPtName;
  
  int step;
  SbVec3f prevPt;
  SbVec3f prevParticleDir;
  void* prevColorField;
  
  viewPtData camB4Animation;
  bool returnToSaveVP;
  bool returnToAnim;
  SoCamera* myCam;
  void setStartingPtForAnimation(); 
  float left_right, up_down;    
  SbVec3f rotAxis; // For 90 degree rotations
  int rotCnt;  // # of steps rotation is split into

  static void getViewPtNameCB(Widget, XtPointer, XtPointer);
  std::string viewPtAutoName();

  ////////////////////////ANIM_SPEED_INDICATOR///////////////////////
  SoSearchAction * searcher;

  SoNode * superimposition;
  SoCoordinate3 * sgeometry;
  SoScale * sscale;
  
  SoTranslation * stranslation;
  SoTranslation * curInfoTrans;
  SoTranslation * mouseOverTransSolid;
  SoTranslation * mouseOverTransMaterial;
  SoTranslation * mouseOverTransLogName;
  SoTranslation * mouseOverTransZPos;

  // Used for 2 similar purposes: 1. Displays z position during animation
  //                              2. Displays name of the current viewpoint
  SoText2 * curInfoText;
  /* Need to use many different fields for mouseover
   * because newlines are ignored when the scene is rendered */
  SoText2 * mouseOverTextSolid;
  SoText2 * mouseOverTextMaterial;
  SoText2 * mouseOverTextLogName;
  SoText2 * mouseOverTextZPos;
   
  SoFont * curInfoFont;
  SoFont * mouseOverFontSolid;
  SoFont * mouseOverFontMaterial;
  SoFont * mouseOverFontLogName;
  SoFont * mouseOverFontZPos;
  SoSwitch * axisSwitch;
  SoSwitch * animSpeedOutlineSwitch; 
  SoSwitch * animSpeedSwitch;
  SoSwitch * curInfoSwitch;
    
  SoNode * getSuperimpositionNode(SoNode *, const char * name);

  void superimpositionEvent(SoAction * action);
  static void superimpositionCB(void * closure, SoAction * action);

  virtual void actualRedraw(void);
  void updateSpeedIndicator(void);

  float maxSpeed; 
  ////////////////////////ANIM_SPEED_INDICATOR///////////////////////

  // FWJ added for Ortho camera
  float defaultHeight;
  float defaultHeightAngle;
  // FWJ add look-ahead for animation tracking on curves
  G4int pathLookahead;
  
  // Used by G4 app during element rotations, stores previous view
  SbVec3f upVector, offsetFromCenter, center;   
  bool rotUpVec;

  SoSeparator* newSceneGraph;


};
#endif /* G4OPENINVENTORXTEXAMINERVIEWER_HH */
