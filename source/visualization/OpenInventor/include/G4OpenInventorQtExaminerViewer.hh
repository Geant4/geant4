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

// Frederick Jones TRIUMF 07 January 2018


#ifndef G4OPENINVENTORQTEXAMINERVIEWER_HH
#define G4OPENINVENTORQTEXAMINERVIEWER_HH

// Set up notification of event processing

#include "G4VStateDependent.hh"

class G4OpenInventorQtExaminerViewer;

class HookEventProcState : public G4VStateDependent
{
public:
   HookEventProcState(G4OpenInventorQtExaminerViewer*);
   ~HookEventProcState();
   virtual G4bool Notify(G4ApplicationState requestedState);
private:
   G4OpenInventorQtExaminerViewer* viewer;
};


#include "G4String.hh"

//#include "G4OpenInventorViewer.hh"

#include <map>
#include <vector>
#include <fstream>
#include <Inventor/SbLinear.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/events/SoKeyboardEvent.h>

#include <qobject.h>

class G4UIQt;

class SoCoordinate3;
class SoFont;
class SoText2;
class SoPointSet;

class QWidget;
class QDialog;
class QMenuBar;
class QMenu;
class QAction;
class QListWidgetItem;
class QPushButton;
class QRadioButton;
class QMessageBox;
class QFont;

// The Aux Window dialog created with Qt Designer:
class Ui_Dialog;


class G4OpenInventorQtExaminerViewer: public QObject,
                                      public SoQtExaminerViewer {

  Q_OBJECT 

  //  friend class G4OpenInventorQtExaminerViewerMessenger;
  // FWJ
 friend class G4OpenInventorQtViewer;

private Q_SLOTS :

   void FileOpenBookmarkCB();
   void FileNewBookmarkCB();
   void FileLoadRefPathCB();
   void FileSaveRefPathCB();
   void FileLoadSceneGraphCB();
   void FileSaveSceneGraphCB();

   void ToolsAnimateRefParticleCB();
   void ToolsRefPathStartCB();
   void ToolsRefPathInvertCB();

   void HelpControlsCB();

   // For added viewer buttons
   void SaveViewPtCB();
   void NextViewPtCB();
   void PrevViewPtCB();
   void AbbrOutputCB(bool);      // Includes mouse-over fcns
   void PickRefPathCB();
   void SwitchWireFrameCB(bool);
   void SwitchAxesCB(bool);
   void DetachCB();

   // Lists Window
   void LoadBookmarkCB(QListWidgetItem*);
   void DeleteBookmarkCB();
   void RenameBookmarkCB();
   void SortBookmarksCB();
   void LookAtSceneElementCB(QListWidgetItem*);

private:

   static G4OpenInventorQtExaminerViewer* viewer;
   QString* fName;

   int OWwidth, OWheight;

   void (*escapeCallback)();

   void* examinerObject;
   SbBool lshiftdown, rshiftdown, lctrldown, rctrldown;

   QFont* font;
   QMenuBar* menubar;
   QMenu* filemenu;
   QMenu* toolsmenu;
   QMenu* etcmenu;
   QMenu* helpmenu;
   QMessageBox* helpmsgbox;

   bool externalQtApp;

   QAction* FileOpenBookmark;
   QAction* FileNewBookmark;
   QAction* FileLoadRefPath;
   QAction* FileSaveRefPath;
   QAction* FileLoadSceneGraph;
   QAction* FileSaveSceneGraph;

   QAction* ToolsAnimateRefParticle;
   QAction* ToolsRefPathStart;
   QAction* ToolsRefPathInvert;

   // KEEP in the viewer
   QAction* HelpControls;

   // Added viewer buttons
   QPushButton* saveViewPtButton;
   QPushButton* nextViewPtButton;
   QPushButton* prevViewPtButton;
   QPushButton* abbrOutputButton;
   QPushButton* pickRefPathButton;
   QPushButton* switchWireFrameButton;
   QPushButton* switchAxesButton;
   QPushButton* detachButton;

   QListWidgetItem* saveViewPtItem;

   Ui_Dialog* AuxWindowDialog;
   QDialog* AuxWindow;

   G4UIQt* uiQt;
   QWidget* viewerParent;
   QWidget* viewerParent2;
   int uiQtTabIndex;

   int processSoEventCount;
   G4String empty = "";

public:

   G4OpenInventorQtExaminerViewer(QWidget* parent = NULL,
             const char* name = NULL,
	     SbBool embed = TRUE, 
	     SoQtFullViewer::BuildFlag flag = BUILD_ALL,
                                  SoQtViewer::Type type = BROWSER);

   ~G4OpenInventorQtExaminerViewer();

   template <class T> void parseString(T &t, const std::string &s, bool &error);

   // In case the viewer is embedded and then detached:
   void setOrigWindowSize(int w, int h) { OWwidth = w; OWheight = h; }

   // Menubar information needed by G4OpenInventorQtViewer
   // for common menu items:
   QMenuBar* getMenubar() { return menubar; }
   QMenu* getFileMenu() { return filemenu; }
   QMenu* getEtcMenu() { return etcmenu; }
   QFont* getFont() { return font; };

   void setExternalQtApp() { externalQtApp = TRUE; }

   // Needed?
   std::string saveScenegraphFileName;
   std::string saveRefCoordsFileName;

   void addEscapeCallback(void (*cb)());

   bool abbrOutputFlag;
   bool pickRefPathFlag;
   bool viewingBeforePickRef;


protected:
   // FWJ Constructor with build flag added (as in parent)
   // Need for this TBD.
   //  G4OpenInventorQtExaminerViewer(QWidget parent,
   //	     const char *name,
   //	     SbBool embed,
   //	     SoQtFullViewer::BuildFlag flag,
   //	     SoQtViewer::Type type,
   //	     SbBool build);

   void construct(const SbBool build);
   void buildWidget(QWidget* parent);

   virtual void afterRealizeHook();

   HookEventProcState* hookBeamOn;
   friend class HookEventProcState;
   bool newEvents;

   static void sceneChangeCB(void*, SoSensor*);

   SbBool processSoEvent(const SoEvent* const event);

   void saveViewPt(char* name);
   bool loadViewPts();
   void addViewPoints();
   void setViewPt();
   void writeViewPtIdx();
   void cleanUpAfterPrevFile();
   void deleteViewPt(char *vpName = NULL);
   void renameViewPt(char *vpName);
   void sortViewPts(std::vector<std::string>);

   void zoom(const float);
   void moveCamera(float dist = 0, bool lookdown = false);
   std::string curEltName;
   SbVec3f camUpVec;
   SbVec3f camDir;
   void rotateCamera();
   void updateViewParams(SoKeyboardEvent::Key);
  
   static void mouseoverCB(void *aThis, SoEventCallback *eventCB);
   static void pickingCB(void *aThis, SoEventCallback *eventCB);


   // Animation
   static void animateRefParticleCB();
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

   void setReferencePath(SoLineSet*, SoCoordinate3*, bool append = false);
   void setReferencePathZPos();
   void findAndSetRefPath();
   SoCoordinate3* getCoordsNode(SoFullPath *path);
   void getSceneElements(); // reads elements from the scene graph
   float sqrlen(const SbVec3f&);
   void distanceToTrajectory(const SbVec3f&, float&, SbVec3f&, int&);
   void sortElements();
   void createElementsList();
   //   static void closeMainWindowCB(Widget, XtPointer, XtPointer);
   void evenOutRefParticlePts();

   //  static void gotoRefPathStartCB(Widget, XtPointer, XtPointer);
   void gotoRefPathStart();
   //  static void invertRefPathCB(Widget, XtPointer, XtPointer);
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

      G4bool operator<(elementForSorting const &other) const
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

   // Need to use many different fields for mouseover
   // because newlines are ignored when the scene is rendered
   SoText2* mouseOverTextSolid;
   SoText2* mouseOverTextMaterial;
   SoText2* mouseOverTextLogName;
   SoText2* mouseOverTextZPos;

   SoFont* curInfoFont;
   SoFont* mouseOverFontSolid;
   SoFont* mouseOverFontMaterial;
   SoFont* mouseOverFontLogName;
   SoFont* mouseOverFontZPos;
   SoSwitch* axisSwitch;
   SoSwitch* animSpeedOutlineSwitch;
   SoSwitch* animSpeedSwitch;
   SoSwitch* curInfoSwitch;

   SoNode* getSuperimpositionNode(SoNode*, const char* name);

   void superimpositionEvent(SoAction* action);
   static void superimpositionCB(void* closure, SoAction* action);


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

#endif /* G4OPENINVENTORQTEXAMINERVIEWER_HH */
