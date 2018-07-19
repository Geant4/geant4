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
// Open Inventor Xt Extended Viewer - 30 Oct 2012
// Rastislav Ondrasek, Pierre-Luc Gagnon, Frederick Jones TRIUMF

#include <stdio.h>
#include <string.h>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <algorithm> // For using sort on a vector
#include <X11/keysym.h>

#include <Xm/Xm.h>
#include <Xm/Text.h>
#include <Xm/List.h>
#include <Xm/MessageB.h>
#include <Xm/PushB.h>
#include <Xm/ToggleB.h>
#include <Xm/CascadeB.h>
#include <Xm/ArrowBG.h>
#include <Xm/Form.h>
#include <Xm/RowColumn.h>
#include <Xm/FileSB.h>
#include <Xm/SelectioB.h>
#include <Xm/Protocols.h>
#include <Xm/SeparatoG.h>
#include <Xm/DialogS.h>  
#include <Xm/PanedW.h>
#include <Xm/LabelG.h>
#include <Xm/Scale.h>
#include <Xm/DrawingA.h>

#include <Inventor/Xt/SoXt.h>
//#include <Inventor/Xt/SoXtInternal.h>
#include <Inventor/Xt/SoXtCursor.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/events/SoLocation2Event.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
// FWJ moved to header file
//#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoWriteAction.h>

#include "G4OpenInventorXtExaminerViewer.hh"
#include "wheelmouse.h"  // To use mouse scrolling in dialogs
#include "SoXtInternal.h"
#include <Inventor/sensors/SoTimerSensor.h>   // Animation
#include <Inventor/sensors/SoNodeSensor.h>    // Detect start of run
#include "saveViewPt.h"
#include "pickext.h"
#include "pickref.h"
#include "wireframe.h"
//#include "console.h"
//#include "favorites.h"

#include "Geant4_SoPolyhedron.h"
//#include "G4RunManager.hh"
//#include "G4Run.hh"
#include "G4TrajectoryPoint.hh"
#include "G4AttHolder.hh"
#include "G4AttCheck.hh"

#include <Inventor/nodes/SoCallback.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>

#include <Inventor/nodes/SoCoordinate3.h>
// For rendering distance during animation:
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoBaseColor.h>

// For searching for nodes within kits:
#include <Inventor/nodekits/SoBaseKit.h>

// FWJ
#include <Inventor/SbVec3f.h>

G4OpenInventorXtExaminerViewer* G4OpenInventorXtExaminerViewer::viewer = 0;

static const char* thisClassName = "G4OpenInventorXtExaminerViewer";
 
#define MIN_SPEED  2.1        // Lower number means faster
#define START_STEP 0.3
#define SPEED_INDICATOR_STEP 0.045
#define MAX_SPEED_INDICATOR  0.81
// Number of steps 90 degree rotation around an element is split into
#define ROT_CNT 6

// Public constructor
G4OpenInventorXtExaminerViewer::G4OpenInventorXtExaminerViewer(Widget parent,
                const char *name, SbBool embed,
		SoXtFullViewer::BuildFlag flag, SoXtViewer::Type type) :
   SoXtExaminerViewer(parent, name, embed, flag, type, FALSE)
{
// Tell GLWidget not to build just yet
   this->constructor(TRUE);
}

// Protected constructor for classes deriving from this viewer.
G4OpenInventorXtExaminerViewer::G4OpenInventorXtExaminerViewer(Widget parent,
                const char *name, SbBool embed,
                SoXtFullViewer::BuildFlag flag, SoXtViewer::Type type,
                SbBool build) :
   SoXtExaminerViewer(parent, name, embed, flag, type, FALSE)
{
   this->constructor(build);
}

// Called by all constructors to set up widgets and initialize member fields.
void G4OpenInventorXtExaminerViewer::constructor(const SbBool build)
{
   setClassName(thisClassName);

   hookBeamOn = new HookEventProcState(this);
   this->newEvents = false;

   fileName = ".bookmarkFile"; // Default viewpoint file name
   viewPtIdx = -1; // index of the most recent viewpoint in viewPtList vector
   animateSensor = new SoTimerSensor(
                       G4OpenInventorXtExaminerViewer::animateSensorCB, this);
   animateSensorRotation = new SoTimerSensor(
               G4OpenInventorXtExaminerViewer::animateSensorRotationCB, this);
   animateBtwPtsPeriod = MIN_SPEED;
   currentState = GENERAL;
   myCam = new SoPerspectiveCamera;
   MAX_VP_IDX = 3;
   MAX_VP_NAME = 35; // Max length of a viewpoint name, padded with spaces
   rotCnt = ROT_CNT; // For 90 degree rotations
   curViewPtName = new char[MAX_VP_NAME + 1];
   left_right = up_down = 0; // For movements around the beam during animation
   speedStep = START_STEP; // For smoother animation speed increase/decrease
   rotUpVec = false; // Used during scene element rotations
   step = 1;	//By default
   // Used for moving along the beam with the
   // mouse instead of rotating the view
   lshiftdown = rshiftdown = false;
   // Used for rotating the view with the camera
   // staying in place
   lctrldown = rctrldown = false;
   // Used to send abbreviated output to the console when
   abbrOutputFlag = false;
   pickRefPathFlag = false;
   prevColorField = NULL;
   warningFlag = false; // We come from the warning dialog
   viewer = this;
   openFileDialog = newFileDialog = listsDialog = (Widget) NULL;
   loadRefCoordsDialog = saveRefCoordsDialog = NULL;
   loadSceneGraphDialog = saveSceneGraphDialog = NULL;
   myElementList = NULL;
   // FWJ default path look-ahead
   pathLookahead = 5;

   newSceneGraph = NULL;
   zcoordSetFlag = false;

   //////////////////////////SUPERIMPOSED SCENE//////////////////////////
   searcher = NULL;
   // Used in animation; progressively scaled for gradual speed change
   maxSpeed = 0.0f;

   static const char * superimposed[] = {
      "#Inventor V2.1 ascii", "",
      "Separator ",
      "{",
      "	MaterialBinding ",
      "	{",
      "   	value OVERALL",
      "	}",
      "  	OrthographicCamera ",
      "	{",
      "    	height 1",
      "		nearDistance 0",
      "    	farDistance 1",
      "	}",
      "  	DEF soxt->callback Callback { }",
      "  	Separator ",
      "	{",
      "   	DEF soxt->translation Translation ",
      "		{",
      "      		translation 0 0 0",
      "	    }",
      "	    DEF soxt->scale Scale ",
      "		{",
      "      		scaleFactor 1 1 1",
      "	    }",
      "		DEF soxt->geometry Coordinate3 ",
      "		{",
      "		    point ",
      "			[",
      "   	 		-0.81	-0.04	0,	-0.81	0		0,",
      "       		-0.81	0.04 	0,	0    	-0.04 	0,",
      "       		0     	0    	0,  0 	    0.04 	0,",
      "       		0.81 	-0.04 	0,  0.81  	0 	    0,",
      "       		0.81  	0.04 	0,",
      "       		0     	0.02 	0,", // idx 9
      "       		0.81  	0.02 	0,  0.81 	-0.02 	0,",
      "       		0    	-0.02 	0,",
      "       		0     	0.01 	0,", // idx 13
      "       		0.4   	0.01 	0,  0.4  	-0.01 	0,",
      "       		0    	-0.01 	0",
      "			]",
      "		}",
      // current speed indicator (outline)
      "    	DEF soxt->animSpeedOutlineSwitch Switch ",
      "		{",
      "      		whichChild -3",
      "      		Material ",
      "			{",
      " 	     	   emissiveColor 0 0 0",
      "		    }",
      "      		IndexedFaceSet ",
      "			{",
      "	        	coordIndex ",
      "				[",
      "          			12, 11, 10, 9, -1",
      " 	        	]",
      "	        }",
      " 		 }",
      // the coordinate system
      "    	DEF soxt->axisSwitch Switch ",
      "		{",
      "  	    	whichChild -3",
      "      		BaseColor ",
      "			{",
      " 	    	    rgb 1 1 1",
      " 	    	}",
      "      		IndexedLineSet ",
      "			{",
      "  	 	    	coordIndex ",
      "				[",
      "          			0, 2, -1,",
      " 			        3, 5, -1,",
      " 			        6, 8, -1,",
      "          			1, 7, -1",
      " 		        ]",
      " 		    }",
      " 		}",
      // current speed indicator
      "    	DEF soxt->animSpeedSwitch Switch ",
      "		{",
      " 		    whichChild -3",
      "      		Material ",
      "			{",
      "       		emissiveColor 0 1 0",
      "       	}",
      "			IndexedFaceSet ",
      "			{",
      "	    		coordIndex ",
      "				[",
      "          			16, 15, 14, 13, -1",
      "  		    	]",
      " 	    	}",
      "    	}",
      "  	}",
      // For displaying either z position (during animation) or current viewpoint name
      "	DEF soxt->curInfoSwitch Switch ",
      "	{",
      "		whichChild -3",
      "    	DEF soxt->curInfoTrans Translation ",
      "		{",
      "      		translation 10 20 30    ",
      "		}",
      "    	DEF soxt->curInfoFont Font ",
      "		{",
      "      		name defaultFont:Bold",
      "      		size 16",
      "  		}",
      "		DEF soxt->curInfoText Text2 ",
      "		{",
      "      		string Hello",
      "	    }",
      "	}",
      // Need to use different fields for mouseover
      // because newlines are ignored when the scene is rendered
      "	Separator ",
      "	{",
      "    	DEF soxt->mouseOverTransLogName Translation ",
      "		{",
      "      		translation 0 0 0    ",
      "		}",
      "    	DEF soxt->mouseOverFontLogName Font ",
      "		{",
      "      		name defaultFont:Bold",
      "      		size 16",
      "  		}",
      "		DEF soxt->mouseOverTextLogName Text2 { } ",
      "	}",
      "	Separator ",
      "	{",
      "    	DEF soxt->mouseOverTransSolid Translation ",
      "		{",
      "      		translation 0 0 0    ",
      "		}",
      "    	DEF soxt->mouseOverFontSolid Font ",
      "		{",
      "      		name defaultFont:Bold",
      "      		size 16",
      "  		}",
      "		DEF soxt->mouseOverTextSolid Text2 { } ",
      "	}",
      "	Separator ",
      "	{",
      "    	DEF soxt->mouseOverTransMaterial Translation ",
      "		{",
      "      		translation 0 0 0    ",
      "		}",
      "    	DEF soxt->mouseOverFontMaterial Font ",
      "		{",
      "      		name defaultFont:Bold",
      "      		size 16",
      "  		}",
      "		DEF soxt->mouseOverTextMaterial Text2 { } ",
      "	}",
      "	Separator ",
      "	{",
      "    	DEF soxt->mouseOverTransZPos Translation ",
      "		{",
      "      		translation 0 0 0    ",
      "		}",
      "    	DEF soxt->mouseOverFontZPos Font ",
      "		{",
      "      		name defaultFont:Bold",
      "      		size 16",
      "  		}",
      "		DEF soxt->mouseOverTextZPos Text2 { } ",
      "	}",
      "}", NULL
   };

   int i, bufsize;
   for (i = bufsize = 0; superimposed[i]; i++)
      bufsize += strlen(superimposed[i]) + 1;
   char * buf = new char[bufsize + 1];
   for (i = bufsize = 0; superimposed[i]; i++) {
      strcpy(buf + bufsize, superimposed[i]);
      bufsize += strlen(superimposed[i]);
      buf[bufsize] = '\n';
      bufsize++;
   }
   SoInput * input = new SoInput;
   input->setBuffer(buf, bufsize);
   SbBool ok = SoDB::read(input, this->superimposition);
   (void)ok;   // FWJ added to avoid compiler warning
   assert(ok);
   delete input;
   delete[] buf;
   this->superimposition->ref();

   this->sscale = (SoScale *) this->getSuperimpositionNode(
                                 this->superimposition, "soxt->scale");
   this->stranslation = (SoTranslation *) this->getSuperimpositionNode(
                                 this->superimposition, "soxt->translation");
   this->sgeometry = (SoCoordinate3 *) this->getSuperimpositionNode(
                                 this->superimposition, "soxt->geometry");
   this->axisSwitch = (SoSwitch *) this->getSuperimpositionNode(
                                 this->superimposition, "soxt->axisSwitch");
   this->animSpeedOutlineSwitch = (SoSwitch *) this->getSuperimpositionNode(
                       this->superimposition, "soxt->animSpeedOutlineSwitch");
   this->animSpeedSwitch = (SoSwitch *) this->getSuperimpositionNode(
                       this->superimposition, "soxt->animSpeedSwitch");
   this->curInfoSwitch = (SoSwitch *) this->getSuperimpositionNode(
                             this->superimposition, "soxt->curInfoSwitch");
   this->curInfoTrans = (SoTranslation *) this->getSuperimpositionNode(
                             this->superimposition, "soxt->curInfoTrans");
   this->curInfoFont = (SoFont *) this->getSuperimpositionNode(
                             this->superimposition, "soxt->curInfoFont");
   this->curInfoText = (SoText2 *) this->getSuperimpositionNode(
                             this->superimposition, "soxt->curInfoText");
   this->mouseOverTransLogName = (SoTranslation*)this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTransLogName");
   this->mouseOverFontLogName = (SoFont *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverFontLogName");
   this->mouseOverTextLogName = (SoText2 *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTextLogName");
   this->mouseOverTransSolid = (SoTranslation *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTransSolid");
   this->mouseOverFontSolid = (SoFont *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverFontSolid");
   this->mouseOverTextSolid = (SoText2 *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTextSolid");
   this->mouseOverTransMaterial = (SoTranslation*)this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTransMaterial");
   this->mouseOverFontMaterial = (SoFont *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverFontMaterial");
   this->mouseOverTextMaterial = (SoText2 *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTextMaterial");
   this->mouseOverTransZPos = (SoTranslation *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTransZPos");
   this->mouseOverFontZPos = (SoFont *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverFontZPos");
   this->mouseOverTextZPos = (SoText2 *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->mouseOverTextZPos");

   SoCallback * cb = (SoCallback *) this->getSuperimpositionNode(
                   this->superimposition, "soxt->callback");
   cb->setCallback(superimpositionCB, this);

   this->addSuperimposition(this->superimposition);
   this->setSuperimpositionEnabled(this->superimposition, FALSE);
   axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
   animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
   animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);

   /////////////////////\SUPERIMPOSED SCENE///////////////////////////////////

   // Build everything else like the parent viewer does
   if (build) {
      Widget w = buildWidget(getParentWidget());
      setBaseWidget(w);

      // Make this window a little bigger because of the extra buttons
      // FWJ but it is already set to 600x600 by vis/open
      //      setSize(SbVec2s(500, 550));
   }

}


// Static function that returns the pointer to G4OpenInventorXtExaminerViewer
// FWJ DISABLED
//G4OpenInventorXtExaminerViewer *G4OpenInventorXtExaminerViewer::getObject()
//{
//   if (!viewer)
//      new G4OpenInventorXtExaminerViewer();
//   return viewer;
//}


// This method locates a named node in the superimposed or original scene.
SoNode *
G4OpenInventorXtExaminerViewer::getSuperimpositionNode(SoNode *root,
                                                     const char * name)
{
   if (!this->searcher)
      this->searcher = new SoSearchAction;
   searcher->reset();
   searcher->setName(SbName(name));
   searcher->setInterest(SoSearchAction::FIRST);
   searcher->setSearchingAll(TRUE);
   searcher->apply(root);
   assert(searcher->getPath());
   return searcher->getPath()->getTail();
}


void G4OpenInventorXtExaminerViewer::superimpositionCB(void * closure,
                                                     SoAction * action)
{
   if (closure)
      ((G4OpenInventorXtExaminerViewer*)closure)->superimpositionEvent(action);
}


// Renders and positions speed indicator and longitudinal
// distance/viewpoint name on the drawing canvas
void G4OpenInventorXtExaminerViewer::superimpositionEvent(SoAction * action)
{

   if (!action->isOfType(SoGLRenderAction::getClassTypeId()))
      return;
   SbViewportRegion vpRegion =
      ((SoGLRenderAction *) action)->getViewportRegion();
   SbVec2s viewportSize = vpRegion.getViewportSizePixels();

   float aspect = float(viewportSize[0]) / float(viewportSize[1]);
   float factorx = 1.0f / float(viewportSize[1]) * 220.0f;
   float factory = factorx;

   if (aspect > 1.0f) {
      this->stranslation->translation.setValue(SbVec3f(0.0f, -0.4f, 0.0f));
   } else {
      this->stranslation->translation.setValue(
                       SbVec3f(0.0f, -0.4f / aspect, 0.0f));
      factorx /= aspect;
      factory /= aspect;
   }
   if (viewportSize[0] > 500)
      factorx *= 500.0f / 400.0f;
   else
      factorx *= float(viewportSize[0]) / 400.0f;
   this->sscale->scaleFactor.setValue(SbVec3f(factorx, factory, 1.0f));

   float xInfo, yInfo, xMouseLogName, yMouseLogName, xMouseSolid, yMouseSolid,
      xMouseMaterial, yMouseMaterial, xMouseZPos, yMouseZPos;
   xInfo = -.45;
   yInfo = .45;
   xMouseLogName = 0.0;
   yMouseLogName = -.75;
   xMouseSolid = 0.0;
   yMouseSolid = -.78;
   xMouseMaterial = 0.0;
   yMouseMaterial = -.81;
   xMouseZPos = 0.0;
   yMouseZPos = -.84;

   if (aspect > 1.0f) {
      xInfo *= aspect;
      xMouseSolid *= aspect;
      xMouseMaterial *= aspect;
      this->curInfoTrans->translation.setValue(SbVec3f(xInfo, yInfo, 0.0));
      this->mouseOverTransLogName->translation.setValue(
                        SbVec3f(xMouseLogName, yMouseLogName, 0.0));
      this->mouseOverTransSolid->translation.setValue(
                        SbVec3f(xMouseSolid, yMouseSolid, 0.0));
      this->mouseOverTransMaterial->translation.setValue(
                        SbVec3f(xMouseMaterial, yMouseMaterial, 0.0));
      this->mouseOverTransZPos->translation.setValue(
                        SbVec3f(xMouseZPos, yMouseZPos, 0.0));
   } else {
      yInfo /= aspect;
      yMouseSolid /= aspect;
      yMouseMaterial /= aspect;
      this->curInfoTrans->translation.setValue(SbVec3f(xInfo, yInfo, 0.0));
      this->mouseOverTransLogName->translation.setValue(
                        SbVec3f(xMouseLogName, yMouseLogName, 0.0));
      this->mouseOverTransSolid->translation.setValue(
                        SbVec3f(xMouseSolid, yMouseSolid, 0.0));
      this->mouseOverTransMaterial->translation.setValue(
                        SbVec3f(xMouseMaterial, yMouseMaterial, 0.0));
      this->mouseOverTransZPos->translation.setValue(
                        SbVec3f(xMouseZPos, yMouseZPos, 0.0));
   }

   if (currentState == VIEWPOINT) { // Displaying viewpoint name
      this->curInfoFont->size.setValue(15);
      this->curInfoFont->name.setValue("defaultFont:Italic");
      this->curInfoText->string.setValue(SbString(curViewPtName));
   }
   else if(currentState == GENERAL) { // Displaying longitudinal distance
      this->curInfoFont->size.setValue(16);
      this->curInfoFont->name.setValue("defaultFont:Bold");
      this->curInfoText->string.setValue(SbString(""));
   }
   else {
      if (refParticleIdx < (int) refParticleTrajectory.size() - 1) {
         this->curInfoFont->size.setValue(16);
         this->curInfoFont->name.setValue("defaultFont:Bold");
         char zPos[20];
         sprintf(zPos, "%7.2f [m]", refZPositions[refParticleIdx] / 1000);
         this->curInfoText->string.setValue(SbString(zPos));
      }
   }
}


G4OpenInventorXtExaminerViewer::~G4OpenInventorXtExaminerViewer()
{
   if (superimposition != NULL) {
      removeSuperimposition(superimposition);
      superimposition->unref();
      superimposition = NULL;
   }
   if (animateSensor->isScheduled())
      animateSensor->unschedule();
   delete animateSensor;
   delete sceneChangeSensor;

   delete[] curViewPtName;
   delete searcher;

   viewer = 0;
}


// Adds a menu bar and a few menu items to the viewer.
Widget G4OpenInventorXtExaminerViewer::buildWidget(Widget parent)
{
   Widget shell;
   Atom WM_DELETE_WINDOW;

   if (!parent)
      SoDebugError::post("G4OpenInventorXtExaminerViewer::buildWidget", "Error: Parent is null.");

   Arg args[10];
   XtSetArg(args[0], XmNtopAttachment, XmATTACH_FORM);
   XtSetArg(args[1], XmNleftAttachment, XmATTACH_FORM);
   XtSetArg(args[2], XmNrightAttachment, XmATTACH_FORM);
   XtSetArg(args[3], XmNbottomAttachment, XmATTACH_FORM);
   Widget form = XmCreateForm(parent, (char *) "Form", args, 4);
   XtManageChild(form);

   shell = XtParent(form);
   WM_DELETE_WINDOW = XInternAtom(XtDisplay(parent), "WM_DELETE_WINDOW",
                                  False);
   XmAddWMProtocolCallback(shell, WM_DELETE_WINDOW,
                           (XtCallbackProc)closeMainWindowCB, this);

   XtSetArg(args[0], XmNtopAttachment, XmATTACH_FORM);
   XtSetArg(args[1], XmNleftAttachment, XmATTACH_FORM);
   XtSetArg(args[2], XmNrightAttachment, XmATTACH_FORM);
   menuBar = XmCreateMenuBar(form, (char *) "MenuBar", args, 3);
   XtManageChild(menuBar);

   fileMenu = addMenu("File");
   this->addButton(fileMenu, "Open Viewpoint File...", openViewPtFileCB);
   addButton(fileMenu, "New Viewpoint File", newViewPtFileCB);
   addButton(fileMenu, "Load Ref. Coords", loadRefCoordsDialogCB);
   addButton(fileMenu, "Save Ref. Coords", saveRefCoordsDialogCB);
   addButton(fileMenu, "Load Scene Graph", loadSceneGraphDialogCB);
   addButton(fileMenu, "Save Scene Graph", saveSceneGraphDialogCB);
   XtManageChild(
          XmCreateSeparatorGadget(fileMenu, (char *) "Separator", NULL, 0));

   Widget menu = addMenu("Tools");
   addButton(menu, "Animate Ref. Particle", animateRefParticleCB);
   addButton(menu, "Go to start of Ref path", gotoRefPathStartCB);
   addButton(menu, "Invert Ref path", invertRefPathCB);

   Widget viewerBase = SoXtFullViewer::buildWidget(form);

   XtSetArg(args[0], XmNtopAttachment, XmATTACH_WIDGET);
   XtSetArg(args[1], XmNtopWidget, menuBar);
   XtSetArg(args[2], XmNleftAttachment, XmATTACH_FORM);
   XtSetArg(args[3], XmNrightAttachment, XmATTACH_FORM);
   XtSetArg(args[4], XmNbottomAttachment, XmATTACH_FORM);
   XtSetValues(viewerBase, args, 5);

   return viewerBase;
}


// Adds a new menu to menuBar
Widget G4OpenInventorXtExaminerViewer::addMenu(std::string name)
{
   Arg args[1];
   Widget menu = XmCreatePulldownMenu(menuBar, (char *) name.c_str(), NULL, 0);

   XtSetArg(args[0], XmNsubMenuId, menu);
   Widget w = XmCreateCascadeButton(menuBar, (char *) name.c_str(), args, 1);
   XtManageChild(w);

   return menu;
}


// Adds a new button to menu
void G4OpenInventorXtExaminerViewer::addButton(Widget menu, std::string name,
                                             XtCallbackProc cb)
{
   Widget button = XmCreatePushButton(menu, (char *) name.c_str(), NULL, 0);
   XtManageChild(button);
   XtAddCallback(button, XmNactivateCallback, cb, this);
}


// Overloaded for saving of and browsing through viewpoints.
void G4OpenInventorXtExaminerViewer::createViewerButtons(Widget parent,
                                                       SbPList * buttonlist)
{
   int n;
   Arg args[6];
   Widget saveViewPtButton, abbrOutputButton, pickRefPathButton;
   Widget switchWireFrameButton;

   // Create original buttons
   SoXtExaminerViewer::createViewerButtons(parent, buttonlist);

   // Handle disappearing button caused by SoXtExaminerViewer::setCamera
   Widget emptyButton = XtVaCreateManagedWidget("", xmPushButtonWidgetClass,
                                                parent, NULL);
   buttonlist->append(emptyButton);

   // Left arrow that goes back one view point on click
   n = 0;
   XtSetArg(args[n], XmNtopPosition, 1);	n++;
   XtSetArg(args[n], XmNbottomPosition, 2);	n++;
   XtSetArg(args[n], XmNleftPosition, 0);	n++;
   XtSetArg(args[n], XmNrightPosition, 1);	n++;
   XtSetArg(args[n], XmNarrowDirection, XmARROW_LEFT);	n++;
   XtSetArg(args[n], XmNsensitive, False);	n++;
   prevViewPtButton = XmCreateArrowButtonGadget(parent, (char *) "ArrowL",
                                                args, n);
   XtManageChild(prevViewPtButton);
   XtAddCallback(prevViewPtButton, XmNactivateCallback,
                 G4OpenInventorXtExaminerViewer::prevViewPtCB, this);
   buttonlist->append(prevViewPtButton);

   // Right arrow that goes forward one view point on click
   n = 0;
   XtSetArg(args[n], XmNtopPosition, 1);	n++;
   XtSetArg(args[n], XmNbottomPosition, 2);	n++;
   XtSetArg(args[n], XmNleftPosition, 0);	n++;
   XtSetArg(args[n], XmNrightPosition, 1);	n++;
   XtSetArg(args[n], XmNarrowDirection, XmARROW_RIGHT);	n++;
   XtSetArg(args[n], XmNsensitive, False);	n++;
   nextViewPtButton = XmCreateArrowButtonGadget(parent, (char *) "ArrowR",
                                                args, n);
   XtManageChild(nextViewPtButton);
   XtAddCallback(nextViewPtButton, XmNactivateCallback,
                 G4OpenInventorXtExaminerViewer::nextViewPtCB, this);
   buttonlist->append(nextViewPtButton);

   // Save button for storing current camera parameters
   saveViewPtButton = XtVaCreateManagedWidget("Save", xmPushButtonWidgetClass,
                                              parent, NULL);
   XtAddCallback(saveViewPtButton, XmNactivateCallback,
                 G4OpenInventorXtExaminerViewer::saveViewPtCB, this);
   Pixmap saveVP, saveVP_ins;
   saveVP = SoXtInternal::createPixmapFromXpm(saveViewPtButton,
                                              saveViewPt_xpm);
   saveVP_ins = SoXtInternal::createPixmapFromXpm(saveViewPtButton,
                                                  saveViewPt_xpm, TRUE);
   XtVaSetValues(saveViewPtButton, XmNlabelType, XmPIXMAP, XmNlabelPixmap,
                 saveVP, XmNselectPixmap, saveVP, XmNlabelInsensitivePixmap,
                 saveVP_ins, XmNselectInsensitivePixmap, saveVP_ins, NULL);
   buttonlist->append(saveViewPtButton);

   // Toggle button to get abbreviated output
   abbrOutputButton = XtVaCreateManagedWidget("Abbr",
                                              xmToggleButtonWidgetClass, parent, XmNindicatorOn, False, NULL);
   XtAddCallback(abbrOutputButton, XmNdisarmCallback, G4OpenInventorXtExaminerViewer::abbrOutputCB,
                 this);
   Pixmap pickextxpm, pickextxpm_ins;
   pickextxpm = SoXtInternal::createPixmapFromXpm(abbrOutputButton,
                                                  pickext_xpm);
   pickextxpm_ins = SoXtInternal::createPixmapFromXpm(abbrOutputButton,
                                                      pickext_xpm, TRUE);
   XtVaSetValues(abbrOutputButton, XmNlabelType, XmPIXMAP, XmNlabelPixmap,
                 pickextxpm, XmNselectPixmap, pickextxpm, XmNlabelInsensitivePixmap,
                 pickextxpm_ins, XmNselectInsensitivePixmap, pickextxpm_ins, NULL);
   //   Pixmap consolexpm, consolexpm_ins;
   // consolexpm = SoXtInternal::createPixmapFromXpm(abbrOutputButton,
   //                                                console_xpm);
   // consolexpm_ins = SoXtInternal::createPixmapFromXpm(abbrOutputButton,
   //                                                    console_xpm, TRUE);
   // XtVaSetValues(abbrOutputButton, XmNlabelType, XmPIXMAP, XmNlabelPixmap,
   //               consolexpm, XmNselectPixmap, consolexpm, XmNlabelInsensitivePixmap,
   //               consolexpm_ins, XmNselectInsensitivePixmap, consolexpm_ins, NULL);
   buttonlist->append(abbrOutputButton);

   // Button for selecting the beam that will act as reference path
   pickRefPathButton = XtVaCreateManagedWidget("Refpath", xmPushButtonWidgetClass,
                                               parent, NULL);
   XtAddCallback(pickRefPathButton, XmNactivateCallback,
                 G4OpenInventorXtExaminerViewer::pickRefPathCB, this);
   Pixmap pickrefxpm, pickrefxpm_ins;
   pickrefxpm = SoXtInternal::createPixmapFromXpm(pickRefPathButton,
                                                    pickref_xpm);
   pickrefxpm_ins = SoXtInternal::createPixmapFromXpm(pickRefPathButton,
                                                        pickref_xpm, TRUE);
   XtVaSetValues(pickRefPathButton, XmNlabelType, XmPIXMAP, XmNlabelPixmap,
       pickrefxpm, XmNselectPixmap, pickrefxpm, XmNlabelInsensitivePixmap,
       pickrefxpm_ins, XmNselectInsensitivePixmap, pickrefxpm_ins, NULL);

   buttonlist->append(pickRefPathButton);

   // Toggle button for switching in and out of wireframe mode
   switchWireFrameButton = XtVaCreateManagedWidget("Wireframe",
         xmToggleButtonWidgetClass, parent,  XmNindicatorOn, False, NULL);
   XtAddCallback(switchWireFrameButton, XmNvalueChangedCallback,
                 G4OpenInventorXtExaminerViewer::switchWireFrameCB, this);
   Pixmap wireframe, wireframe_ins;
   wireframe = SoXtInternal::createPixmapFromXpm(switchWireFrameButton,
                                                 wireframe_xpm);
   wireframe_ins = SoXtInternal::createPixmapFromXpm(switchWireFrameButton,
                                                     wireframe_xpm, TRUE);
   XtVaSetValues(switchWireFrameButton, XmNlabelType, XmPIXMAP, XmNlabelPixmap,
              wireframe, XmNselectPixmap, wireframe, XmNlabelInsensitivePixmap,
              wireframe_ins, XmNselectInsensitivePixmap, wireframe_ins, NULL);
   buttonlist->append(switchWireFrameButton);
}


// Called right after buttons and widgets get realized.
// It sets the viewpoint last accessed.
void G4OpenInventorXtExaminerViewer::afterRealizeHook()
{
   SoXtExaminerViewer::afterRealizeHook();

   // Default height is used when selecting and viewing scene elements
   // FWJ Added defaultHeight for Ortho camera
   SoCamera *cam = getCamera();
   if (cam) {
      if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
         defaultHeightAngle =
            ((SoPerspectiveCamera *) cam)->heightAngle.getValue();
         toggleCameraType();
         defaultHeight =
            ((SoOrthographicCamera *) cam)->height.getValue();
         toggleCameraType();
      } else {
         defaultHeight =
            ((SoOrthographicCamera *) cam)->height.getValue();
         toggleCameraType();
         cam = getCamera();
         if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))
            defaultHeightAngle =
               ((SoPerspectiveCamera *) cam)->heightAngle.getValue();
         toggleCameraType();
      }
   }

   // Open the default bookmark file
   fileIn.open(fileName.c_str());
   if (!fileIn.fail()) {
      if (!loadViewPts()) {
         String dialogName = (char *) "Error Loading File";
         std::string msg = "Wrong or corrupted input file.";
         warningMsgDialog(msg, dialogName, NULL);
      } else {
         // Opens a file without erasing it
         fileOut.open(fileName.c_str(), std::ios::in);
         fileOut.seekp(0, std::ios::end); // For appending new data to the end
         constructListsDialog(getParentWidget(), this, NULL); // Pop up listsDialog
         if (viewPtList.size()) {
            // FWJ disabled auto-selection of first viewpoint.
            // Initial view should be user-controllable & not forced
            //    setViewPt();
            XtSetSensitive(nextViewPtButton, True);
            XtSetSensitive(prevViewPtButton, True);
         }
      }
      fileIn.close();
   } else {
      // Creates a new default bookmark file
      fileOut.open(fileName.c_str());
      constructListsDialog(getParentWidget(), this, NULL); // Pop up listsDialog
   }

   fileIn.clear();

   SoSeparator *root = (SoSeparator *) (getSceneManager()->getSceneGraph());
   if (root == NULL)
      SoDebugError::post("G4OpenInventorXtExaminerViewer::afterRealizeHook", "Root is null.");
   else {
      root->addChild(myCam); // For position/orientation calculation during animation
   }

   sceneChangeSensor = new SoNodeSensor;
   sceneChangeSensor->setFunction(sceneChangeCB);
   sceneChangeSensor->attach(root);
   sceneChangeSensor->setData(this);

   // Monitor mouseover events for displaying the name of scene elements
   // An SoEventCallback is needed instead of using the default processSoEvent
   // because that last one does not provide us with an SoPath to the object
   // that was picked
   SoEventCallback *moCB = new SoEventCallback;
   moCB->addEventCallback(
                          SoLocation2Event::getClassTypeId(),
                          mouseoverCB, static_cast<void *>(this));
   root->addChild(moCB);

   // Override the default picking mechanism present in G4OpenInventorViewer
   // because we want abbreviated output when picking a trajectory
   SoEventCallback *pickCB = new SoEventCallback;
   pickCB->addEventCallback(
                            SoMouseButtonEvent::getClassTypeId(),
                            pickingCB, static_cast<void *>(this));
   root->addChild(pickCB);

}


// Rotates camera 90 degrees around a scene element.
// Rotation is animated for smoothness.
void G4OpenInventorXtExaminerViewer::rotateCamera()
{
   SoCamera *cam = getCamera();

   SbRotation rot(rotAxis, M_PI / (2 * ROT_CNT));
   rot.multVec(camDir, camDir);
   rot.multVec(camUpVec, camUpVec);

   SbVec3f camPosNew = prevPt - (camDir*distance);
   cam->position = camPosNew;
   cam->pointAt(prevPt, camUpVec);
   cam->focalDistance = (prevPt - camPosNew).length();

   rotCnt--;

   if (animateSensorRotation->isScheduled()) {
      animateSensorRotation->unschedule();
   }

   animateSensorRotation->setBaseTime(SbTime::getTimeOfDay());
   animateSensorRotation->setInterval(SbTime(0.02));
   animateSensorRotation->schedule();

}


// Slides camera along the beamline.
void G4OpenInventorXtExaminerViewer::moveCamera(float dist, bool lookdown)
{

   SoCamera *cam = getCamera();
   SbVec3f p1, p2;	 // The particle moves from p1 to p2
   SbVec3f particleDir;	 // Direction vector from p1 to p2
   SbVec3f camPosNew;	 // New position of the camera

   if(refParticleTrajectory.size() == 0) {
      //refParticleTrajectory hasn't been set yet
      if(dist)
         distance = dist;
      else
         distance = (cam->position.getValue() - center).length();

      cam->position.setValue(center + offsetFromCenter*distance);
      cam->focalDistance = (cam->position.getValue() - center).length();
      cam->pointAt(center, upVector);
   }
   else{

      // If we move forward past the last trajectory point,
      // go back to the beginning
      if (refParticleIdx >= (int) refParticleTrajectory.size() - 1) {
         prevPt = refParticleTrajectory[refParticleIdx - step];
         dist = (prevPt - cam->position.getValue()).length();
         refParticleIdx = 0;
      }
      // If we move backward past the beginning,
      // go to the last trajectory point
      if (refParticleIdx < 0) {
         prevPt = refParticleTrajectory[refParticleIdx + step];
         dist = (prevPt - cam->position.getValue()).length();
         refParticleIdx = refParticleTrajectory.size() - 2;
      }

      // Set start and end points
      p1 = refParticleTrajectory[refParticleIdx];
      p2 = refParticleTrajectory[refParticleIdx + step];

      // Get the direction from p1 to p2
      particleDir = p2 - p1;
      particleDir.normalize();

      if(prevParticleDir == SbVec3f(0,0,0)){
         // First time entering BEAMLINE mode, look at
         // the element from the front, with camera upright
         if(lookdown)
            camDir = SbVec3f(0,0,1);
         else
            camDir = SbVec3f(1,0,0);
         camUpVec = SbVec3f(0,1,0);

         // In case the start of the goes in a
         // direction other than +z, rotate the camera accordingly
         SbRotation rot(SbVec3f(0,0,1), particleDir);
         rot.multVec(camDir, camDir);
         rot.multVec(camUpVec, camUpVec);

      }
      else if(particleDir != prevParticleDir) {
         // The beamline has changed direction

         SbRotation rot(prevParticleDir, particleDir);
         rot.multVec(camDir, camDir);
         rot.multVec(camUpVec, camUpVec);

      }

      if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
         if (!dist)
            distance = (prevPt - cam->position.getValue()).length();
         else
            distance = dist;
      }

      // FWJ distance not relevant -- use focalDistance
      // if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
      //    if (!dist)
      //       distance = (prevPt - cam->position.getValue()).length();
      //    else
      //       distance = dist;
      // }


      float x,y,z;
      prevPt.getValue(x,y,z);


      if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
         camPosNew = p2 - (camDir*distance);
      }
      if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
         // FWJ maintain focal distance
         camPosNew = p2 - (camDir*cam->focalDistance.getValue());
         //         camPosNew = p2 - (camDir);
      }

      cam->position = camPosNew;
      cam->pointAt(p2, camUpVec);
      cam->focalDistance = (p2 - camPosNew).length();

      p2.getValue(x,y,z);
      camPosNew.getValue(x,y,z);

      prevParticleDir = particleDir;
      prevPt = p1; // For accurate distance calculation

   }

}


void G4OpenInventorXtExaminerViewer::pickingCB(void *aThis, 
                                               SoEventCallback *eventCB)
{
   SoHandleEventAction* action = eventCB->getAction();
   const SoPickedPoint *pp = action->getPickedPoint();
   G4OpenInventorXtExaminerViewer* This = (G4OpenInventorXtExaminerViewer*)aThis;

   if(pp != NULL) {

      SoPath* path = pp->getPath();
      SoNode* node = ((SoFullPath*)path)->getTail();

      if(node->getTypeId() == SoLineSet::getClassTypeId()){

         if(This->pickRefPathFlag){
            This->pickRefPathFlag = false;
            if(This->viewingBeforePickRef != This->isViewing())
               This->setViewing(This->viewingBeforePickRef);
            else
               This->setComponentCursor(SoXtCursor(SoXtCursor::DEFAULT));

            // The trajectory is a set of lines stored in a LineSet
            SoLineSet * trajectory = (SoLineSet *)node;

       // The set of all trajectories is stored in a Seperator group node
       // one level above the LineSet that was picked. The nodes under that
       // seperator are as follows (in this order): Material, LightModel,
       // ResetTransform, MatrixTransform, Coordinate3, DrawStyle, LineSet
            SoSeparator * grpNode = 
               (SoSeparator*)(((SoFullPath*)path)->getNodeFromTail(1));

   // The node that contains the coordinates for the trajectory is a
   // Coordinate3 node which occurs before the LineSet node.  We iterate
   // back through the nodes in the group until we find the Coordinate3 node
            int nodeIndex = grpNode->findChild(trajectory);
            SoNode * tmpNode;
            // FWJ needs initialization
            SoCoordinate3 * coords = 0;
            //            SoCoordinate3 * coords;
            // We allow only 100 iterations, in case the node isn't found
            // (should take only a few iterations)
            for(int i = 0; i < 100; ++i) {
               --nodeIndex;

               tmpNode = grpNode->getChild(nodeIndex);
               if(tmpNode->getTypeId() == SoCoordinate3::getClassTypeId()){
                  //node found
                  coords = (SoCoordinate3 *)tmpNode;
                  break;
               }
            }

            if(coords == NULL){
               String dialogName = (char *) "No coordinates";
               std::string msg = "Could not find the coordinates node"
                  " for the picked trajectory."
                  " Reference trajectory not set";
               This->warningMsgDialog(msg, dialogName, NULL);
               return;
            }


            if ((This->lshiftdown)	|| (This->rshiftdown))
               This->setReferencePath(trajectory, coords, true);
            else
               This->setReferencePath(trajectory, coords, false);

            return;

         }
         else if(This->abbrOutputFlag) {

            G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
            if(attHolder && attHolder->GetAttDefs().size()) {

               std::string strTrajPoint = "G4TrajectoryPoint:";
               std::ostringstream oss;
               for (size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
                  G4cout << G4AttCheck(attHolder->GetAttValues()[i],
                                       attHolder->GetAttDefs()[i]);
                  oss << G4AttCheck(attHolder->GetAttValues()[i],
                                    attHolder->GetAttDefs()[i]);
                  if(oss.str().find(strTrajPoint) != std::string::npos) {

           // Last attribute displayed was a trajectory point.  Since we
           // want abbreviated output, display the last one and exit
           // (unless we're already at the last (and only) trajectory point)
                     if(i != attHolder->GetAttDefs().size()-1) {
                        G4cout << G4AttCheck(
              attHolder->GetAttValues()[attHolder->GetAttDefs().size()-1],
              attHolder->GetAttDefs()[attHolder->GetAttDefs().size()-1]);
                     }
                     break;
                  }
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

            return;
         }
         else{
            //Go to default behavior
         }
      }
      else {
         //Go to default behavior
      }

      // Default behavior in G4OpenInventorViewer::SelectionCB
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

      //Suppress other event handlers
      eventCB->setHandled();
   }
}


void G4OpenInventorXtExaminerViewer::mouseoverCB(void *aThis, SoEventCallback *eventCB)
{
   SoHandleEventAction* action = eventCB->getAction();
   const SoPickedPoint* pp = action->getPickedPoint();
   G4OpenInventorXtExaminerViewer* This = (G4OpenInventorXtExaminerViewer*)aThis;

   if(!This->abbrOutputFlag)
      return;

   if(pp != NULL) {

      const SbViewportRegion & viewportRegion = action->getViewportRegion();

      std::string sLogName;
      float x,y,z;
      std::stringstream ssZPos;
      std::stringstream ssSolids;
      std::stringstream ssMaterials;
      SoPath * path = pp->getPath();
      SoNode* node = ((SoFullPath*)path)->getTail();

      if(node->getTypeId() == Geant4_SoPolyhedron::getClassTypeId()) {

         sLogName = "Logical Volume:  ";
         sLogName += ((Geant4_SoPolyhedron *)node)->getName().getString();

         SoGetBoundingBoxAction bAction(viewportRegion);
         bAction.apply((SoFullPath*)path);
         SbBox3f bBox = bAction.getBoundingBox();
         SbVec3f center = bBox.getCenter();
         center.getValue(x,y,z);
         ssZPos << "Pos:  " << x << "  " << y << "  " << z;

         G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
         if(attHolder && attHolder->GetAttDefs().size()) {

            std::vector<const std::map<G4String,G4AttDef>*> vecDefs =
               attHolder->GetAttDefs();
            std::vector<const std::vector<G4AttValue>*> vecVals =
               attHolder->GetAttValues();
            for (size_t i = 0; i < vecDefs.size(); ++i) {
               const std::vector<G4AttValue> * vals = vecVals[i];

               std::vector<G4AttValue>::const_iterator iValue;

               for (iValue = vals->begin(); iValue != vals->end(); ++iValue) {
                  const G4String& valueName = iValue->GetName();
                  const G4String& value = iValue->GetValue();

                  if(valueName == "Solid") {
                     if(ssSolids.str() == "")
                        ssSolids << "Solid Name:  " << value;
                     else
                        ssSolids << ", " << value;
                  }

                  if(valueName == "Material") {
                     if(ssMaterials.str() == "")
                        ssMaterials << "Material Name:  " << value;
                     else
                        ssMaterials << ", " << value;
                  }
               }
            }
         }
      }
      // FWJ Mouseover for trajectories
      else if(node->getTypeId() == SoLineSet::getClassTypeId()) {
         //         G4cout << "Trajectory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
         G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
         if(attHolder && attHolder->GetAttDefs().size()) {
            std::string strTrajPoint = "G4TrajectoryPoint:";
            std::ostringstream oss;
            G4String t1, t1Ch, t2, t3, t4;
            for (size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
               //               G4cout << "Getting index " << i << " from attHolder" << G4endl;
               // No, returns a vector!   G4AttValue* attValue = attHolder->GetAttValues()[i];
               const std::vector<G4AttValue>* vals = attHolder->GetAttValues()[i];
               std::vector<G4AttValue>::const_iterator iValue;
               for (iValue = vals->begin(); iValue != vals->end(); ++iValue) {
                  const G4String& valueName = iValue->GetName();
                  const G4String& value = iValue->GetValue();
                  // G4cout << "  valueName = " << valueName << G4endl;
                  // G4cout << "  value = " << value << G4endl;
                  // LINE 1
                  if (valueName == "PN") t1 = value;
                  if (valueName == "Ch") {
                     if (atof(value.c_str()) > 0)
                        t1Ch = "    +";
                     else
                        t1Ch = "    ";
                     t1Ch += value;
                  }
                  if (valueName == "PDG") {
                     t1 += "    ";
                     t1 += value;
                     t1 += t1Ch;
                     This->mouseOverTextLogName->string.setValue(t1);
                  }
                  //                  G4cout << "  t1 = " << t1 << G4endl;
                  // LINE 2
                  if (valueName == "EventID") t2 = "Evt " + value;
                  if (valueName == "ID") t2 += "    Trk " + value;
                  if (valueName == "PID") {
                     t2 += "    Prt " + value;
                     This->mouseOverTextSolid->string.setValue(t2);
                  }
                  // LINE 3
                  if (valueName == "IKE") t3 = "KE " + value;
                  if (valueName == "IMom") {
                     // Remove units
                     unsigned ipos = value.rfind(" ");
                     G4String value1 = value;
                     value1.erase(ipos);
                     t3 += "    P (" + value1 + ")";
                  }
                  if (valueName == "IMag") {
                     t3 += " " + value + "/c";
                     //                     t3 += " " + value;
                     This->mouseOverTextMaterial->string.setValue(t3);
                  }
                  // LINE 4
                  if (valueName == "NTP") {
                     std::ostringstream t4oss;
                     t4oss << "TrjPts " <<  value;
                     t4oss << "    Pos " << pp->getPoint()[0] << " " << pp->getPoint()[1] <<
                        " " << pp->getPoint()[2];
                     This->mouseOverTextZPos->string.setValue(SbString(t4oss.str().c_str()));
                  }
               }
//               G4cout << "  NOW CALLING G4AttCheck" << G4endl;
//                G4cout << G4AttCheck(attHolder->GetAttValues()[i],
//                                     attHolder->GetAttDefs()[i]);
//                oss << G4AttCheck(attHolder->GetAttValues()[i],
//                                  attHolder->GetAttDefs()[i]);
//                if(oss.str().find(strTrajPoint) != std::string::npos) {
//                   // Last attribute displayed was a trajectory point.  Since we
//                   // want abbreviated output, display the last one and exit
//                   // (unless we're already at the last (and only) trajectory point)
//                   if(i != attHolder->GetAttDefs().size()-1) {
//                      G4cout << G4AttCheck(
//                                           attHolder->GetAttValues()[attHolder->GetAttDefs().size()-1],
//                                           attHolder->GetAttDefs()[attHolder->GetAttDefs().size()-1]);
//                   }
//                   break;
//                }
            }
         }
         This->setSuperimpositionEnabled(This->superimposition, TRUE);
         This->scheduleRedraw();
         eventCB->setHandled();
         return;
      }

      bool redraw = false;
      if(std::string(This->mouseOverTextLogName->string.getValues(0)->getString()) != sLogName) {
         This->mouseOverTextLogName->string.setValue(SbString(sLogName.c_str()));
         redraw = true;
      }
      if(std::string(This->mouseOverTextSolid->string.getValues(0)->getString()) != ssSolids.str()) {
         This->mouseOverTextSolid->string.setValue(SbString(ssSolids.str().c_str()));
         redraw = true;
      }
      if(std::string(This->mouseOverTextMaterial->string.getValues(0)->getString()) != ssMaterials.str()){
         This->mouseOverTextMaterial->string.setValue(SbString(ssMaterials.str().c_str()));
         redraw = true;
      }
      if(std::string(This->mouseOverTextZPos->string.getValues(0)->getString()) != ssZPos.str()) {
         This->mouseOverTextZPos->string.setValue(SbString(ssZPos.str().c_str()));
         redraw = true;
      }

      if(redraw) {
         This->setSuperimpositionEnabled(This->superimposition, TRUE);
         This->scheduleRedraw();
      }

      eventCB->setHandled();
   }
   else {
      if(std::string(This->mouseOverTextLogName->string.getValues(0)->getString()) != "") {
         This->mouseOverTextLogName->string.setValue(SbString(""));
         This->scheduleRedraw();
      }
      if(std::string(This->mouseOverTextSolid->string.getValues(0)->getString()) != "") {
         This->mouseOverTextSolid->string.setValue(SbString(""));
         This->scheduleRedraw();
      }
      if(std::string(This->mouseOverTextMaterial->string.getValues(0)->getString()) != "") {
         This->mouseOverTextMaterial->string.setValue(SbString(""));
         This->scheduleRedraw();
      }
      if(std::string(This->mouseOverTextZPos->string.getValues(0)->getString()) != "") {
         This->mouseOverTextZPos->string.setValue(SbString(""));
         This->scheduleRedraw();
      }
   }
}


SbBool G4OpenInventorXtExaminerViewer::processSoEvent(const SoEvent * const ev) {
   SoCamera *cam = getCamera();
   const SoType type(ev->getTypeId());

   if (type.isDerivedFrom(SoMouseButtonEvent::getClassTypeId())) {
      SoMouseButtonEvent * me = (SoMouseButtonEvent *) ev;

      if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
          || currentState == PAUSED_ANIMATION) {
         switch (me->getButton()) {
         case SoMouseButtonEvent::BUTTON4: // Scroll wheel up
            if (me->getState() == SoButtonEvent::DOWN) {
               if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
                  float hAngle =
                     ((SoPerspectiveCamera *) cam)->heightAngle.getValue();
                  ((SoPerspectiveCamera *) cam)->heightAngle = hAngle
                     + 0.01f;
                  return TRUE;
               } else if (cam->isOfType(
                                        SoOrthographicCamera::getClassTypeId())) {
                  float height =
                     ((SoOrthographicCamera *) cam)->height.getValue();
                  ((SoOrthographicCamera *) cam)->height = height + 5;
                  return TRUE;
               }
            }
            break;
         case SoMouseButtonEvent::BUTTON5: // Scroll wheel down
            if (me->getState() == SoButtonEvent::DOWN) {
               if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
                  float hAngle =
                     ((SoPerspectiveCamera *) cam)->heightAngle.getValue();
                  if (hAngle > 0.01)
                     ((SoPerspectiveCamera *) cam)->heightAngle = hAngle
                        - 0.01f;
                  return TRUE;
               } else if (cam->isOfType(
                                        SoOrthographicCamera::getClassTypeId())) {
                  float height =
                     ((SoOrthographicCamera *) cam)->height.getValue();
                  if (height > 5)
                     ((SoOrthographicCamera *) cam)->height = height - 5;
                  return TRUE;
               }
            }
            break;
         default:
            break;
         }
      }
      if (currentState == GENERAL) {

      }
   }

   if (type.isDerivedFrom(SoKeyboardEvent::getClassTypeId())) {
      SoKeyboardEvent * ke = (SoKeyboardEvent *) ev;

      if (SoKeyboardEvent::isKeyPressEvent(ev, ke->getKey())) {
         switch (ke->getKey()) {
         case SoKeyboardEvent::LEFT_SHIFT:
            this->lshiftdown = true;
            return TRUE;
         case SoKeyboardEvent::RIGHT_SHIFT:
            this->rshiftdown = true;
            return TRUE;
         case SoKeyboardEvent::LEFT_CONTROL:
            this->lctrldown = true;
            return TRUE;
         case SoKeyboardEvent::RIGHT_CONTROL:
            this->rctrldown = true;
            return TRUE;
         case SoKeyboardEvent::SPACE:
            if (currentState == ANIMATION
                || currentState == REVERSED_ANIMATION) {
               beforePausing = currentState;
               currentState = PAUSED_ANIMATION;
               if (animateSensor->isScheduled())
                  animateSensor->unschedule();
               return TRUE;
            } else if (currentState == PAUSED_ANIMATION) {
               if (maxSpeed) {
                  if ((beforePausing == ANIMATION
                       && refParticleIdx
                       < (int) refParticleTrajectory.size() - 1)
                      || (beforePausing == REVERSED_ANIMATION
                          && refParticleIdx > 0)) {
                     currentState = beforePausing;
                     animateRefParticle();
                  }
               }
               return TRUE;
            }
            break;
         case SoKeyboardEvent::ESCAPE:
            if (currentState == ANIMATION
                || currentState == REVERSED_ANIMATION
                || currentState == PAUSED_ANIMATION) {

               if (animateSensor->isScheduled())
                  animateSensor->unschedule();
               currentState = prevState;
               refParticleIdx = prevRefIdx;
               setSuperimpositionEnabled(superimposition, FALSE);
               maxSpeed = 0.0f;
               step = 1;

               scheduleRedraw();
               if (currentState == VIEWPOINT) {
                  setSuperimpositionEnabled(superimposition, TRUE);
                  axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
                  animSpeedOutlineSwitch->whichChild.setValue(
                                                              SO_SWITCH_NONE);
                  animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);

                  scheduleRedraw();
               }
               restoreCamera();
               return TRUE;
            }
            break;
         case SoKeyboardEvent::DELETE:
            if (viewPtList.size()
                && (currentState != ANIMATION
                    || currentState != REVERSED_ANIMATION
                    || currentState != PAUSED_ANIMATION)) {
               String dialogName = (char *) "Delete Viewpoint";
               std::string msg = "Are you sure you want to delete current viewpoint?";
               warningMsgDialog(msg, dialogName, deleteViewPtCB);
               return TRUE;
            }
            break;
         case SoKeyboardEvent::LEFT_ARROW:
            switch (currentState) {
            case BEAMLINE:
               if ((this->lshiftdown)	|| (this->rshiftdown)){
                  refParticleIdx -= step;
                  moveCamera();
               }
               else if ((this->lctrldown)	|| (this->rctrldown)){
                  if (SoXtExaminerViewer::isAnimating())
                     stopAnimating();
                  prevState = currentState;
                  currentState = ROTATING;
                  animateBtwPtsPeriod = 0.08f;

                  SbVec3f tmp = camDir;
                  tmp.negate();
                  rotAxis = tmp;

                  rotCnt = ROT_CNT;
                  moveCamera(); // To make sure camera is perpendicular to the beamline

                  rotateCamera();
               }
               else{
                  if (SoXtExaminerViewer::isAnimating())
                     stopAnimating();
                  prevState = currentState;
                  currentState = ROTATING;
                  animateBtwPtsPeriod = 0.08f;

                  SbVec3f tmp = camUpVec;
                  tmp.negate();
                  rotAxis = tmp;

                  rotCnt = ROT_CNT;
                  moveCamera(); // To make sure camera is perpendicular to the beamline

                  rotateCamera();

               }
               return TRUE;

            case ANIMATION:
            case REVERSED_ANIMATION:
               left_right -= 1.5f;
               return TRUE;
            case PAUSED_ANIMATION:
               left_right -= 1.5f;
               setStartingPtForAnimation();
               cam->position = myCam->position;
               return TRUE;
            case GENERAL:
            case VIEWPOINT:
               if ((!this->lshiftdown) && (!this->rshiftdown)) {
                  // Using this allows us to look around without
                  // changing the camera parameters (camDir, camUpVec)
                  this->bottomWheelMotion(
                                          this->getBottomWheelValue() + 0.1f);

                  return TRUE;
               }
               break;
            case ROTATING:
               // For this state, let the keyboard event
               // be handled by superclass
               break;
            default:
               SoDebugError::post("G4OpenInventorXtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::RIGHT_ARROW:
            switch(currentState){
            case BEAMLINE:
               if ((this->lshiftdown)	|| (this->rshiftdown)){
                  refParticleIdx += step;
                  moveCamera();
               }
               else if ((this->lctrldown)	|| (this->rctrldown)){
                  if (SoXtExaminerViewer::isAnimating())
                     stopAnimating();
                  prevState = currentState;
                  currentState = ROTATING;
                  animateBtwPtsPeriod = 0.08f;

                  rotAxis = camDir;

                  rotCnt = ROT_CNT;
                  moveCamera(); // To make sure camera is perpendicular to the beamline

                  rotateCamera();
               }
               else{
                  if (SoXtExaminerViewer::isAnimating())
                     stopAnimating();
                  prevState = currentState;
                  currentState = ROTATING;
                  animateBtwPtsPeriod = 0.08f;

                  rotAxis = camUpVec;

                  rotCnt = ROT_CNT;
                  moveCamera(); // To make sure camera is perpendicular to the beamline

                  rotateCamera();
               }
               return TRUE;

            case ANIMATION:
            case REVERSED_ANIMATION:
               left_right += 1.5f;
               return TRUE;
            case PAUSED_ANIMATION:
               left_right += 1.5f;
               setStartingPtForAnimation();
               cam->position = myCam->position;
               return TRUE;
            case GENERAL:
            case VIEWPOINT:
               if ((!this->lshiftdown) && (!this->rshiftdown)) {
                  // Using this allows us to look around without
                  // changing the camera parameters (camDir, camUpVec)
                  this->bottomWheelMotion(
                                          this->getBottomWheelValue() - 0.1f);
                  return TRUE;
               }
               break;
            case ROTATING:
               // For this state, let the keyboard event
               // be handled by superclass
               break;
            default:
               SoDebugError::post("G4OpenInventorXtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::DOWN_ARROW:
            switch(currentState){
            case BEAMLINE:

               if ((this->lshiftdown)	|| (this->rshiftdown)){
                  refParticleIdx -= step;
                  moveCamera();
               }
               else{
                  if (SoXtExaminerViewer::isAnimating())
                     stopAnimating();
                  prevState = currentState;
                  currentState = ROTATING;
                  animateBtwPtsPeriod = 0.08f;

                  rotAxis = camDir.cross(camUpVec);

                  rotCnt = ROT_CNT;
                  moveCamera(); // To make sure camera is perpendicular to the beamline

                  rotateCamera();

               }
               return TRUE;

            case ANIMATION:
            case REVERSED_ANIMATION:
               up_down -= 1.5f;
               return TRUE;
            case PAUSED_ANIMATION:
               up_down -= 1.5f;
               setStartingPtForAnimation();
               cam->position = myCam->position;
               return TRUE;
            case GENERAL:
            case VIEWPOINT:
               // Using this allows us to look around without
               // changing the camera parameters (camDir, camUpVec)
               if ((!this->lshiftdown) && (!this->rshiftdown)) {
                  this->leftWheelMotion(this->getLeftWheelValue() - 0.1f);
                  return TRUE;
               }
               break;
            case ROTATING:
               // For this state, let the keyboard event
               // be handled by superclass
               break;
            default:
               SoDebugError::post("G4OpenInventorXtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::UP_ARROW:
            switch(currentState){
            case BEAMLINE:
               if ((this->lshiftdown)	|| (this->rshiftdown)){
                  refParticleIdx -= step;
                  moveCamera();
               }
               else{
                  if (SoXtExaminerViewer::isAnimating())
                     stopAnimating();
                  prevState = currentState;
                  currentState = ROTATING;
                  animateBtwPtsPeriod = 0.08f;

                  rotAxis = camUpVec.cross(camDir);

                  rotCnt = ROT_CNT;
                  moveCamera();

                  rotateCamera();


               }
               return TRUE;
            case ANIMATION:
            case REVERSED_ANIMATION:
               up_down += 1.5f;
               return TRUE;
            case PAUSED_ANIMATION:
               up_down += 1.5f;
               setStartingPtForAnimation();
               cam->position = myCam->position;
               return TRUE;
            case GENERAL:
            case VIEWPOINT:
               // Using this allows us to look around without
               // changing the camera parameters (camDir, camUpVec)
               if ((!this->lshiftdown) && (!this->rshiftdown)) {
                  this->leftWheelMotion(this->getLeftWheelValue() + 0.1f);
                  return TRUE;
               }
               break;
            case ROTATING:
               // For this state, let the keyboard event
               // be handled by superclass
               break;
            default:
               SoDebugError::post("G4OpenInventorXtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::PAGE_UP:
            switch(currentState){
            case BEAMLINE:
               if (step < (int) refParticleTrajectory.size() / 5) // Magic number
                  step++;
               return TRUE;
            case ANIMATION:
               incSpeed();
               maxSpeed += SPEED_INDICATOR_STEP;
               if (maxSpeed > 0.8)
                  maxSpeed = MAX_SPEED_INDICATOR;
               scheduleRedraw();

               return TRUE;
            case REVERSED_ANIMATION:
               if(!animateSensor->isScheduled()){
                  currentState = ANIMATION;
                  if (refParticleIdx
                      < (int) refParticleTrajectory.size() - 1) {
                     refParticleIdx++;
                     maxSpeed = SPEED_INDICATOR_STEP;
                     scheduleRedraw();
                     animateRefParticle();
                  }
               }
               else{
                  maxSpeed += SPEED_INDICATOR_STEP;
                  decSpeed();
                  scheduleRedraw();
               }
               return TRUE;
            case PAUSED_ANIMATION:
               maxSpeed += SPEED_INDICATOR_STEP;
               if (maxSpeed > 0.8)
                  maxSpeed = MAX_SPEED_INDICATOR;

               if (beforePausing == ANIMATION) {
                  incSpeed();
               } else {
                  decSpeed();
                  if (animateBtwPtsPeriod >= MIN_SPEED)
                     beforePausing = ANIMATION;
               }

               scheduleRedraw();
               return TRUE;
            default:	//fall through
               break;
            }
            break;

         case SoKeyboardEvent::PAGE_DOWN:
            switch(currentState){
            case BEAMLINE:
               if (step > 1)
                  step--;
               return TRUE;
            case ANIMATION:
               if(!animateSensor->isScheduled()){
                  currentState = REVERSED_ANIMATION;
                  if (refParticleIdx > 1) {
                     refParticleIdx--;
                     maxSpeed = -SPEED_INDICATOR_STEP;
                     scheduleRedraw();
                     animateRefParticle();
                  }
               }
               else{
                  maxSpeed -= SPEED_INDICATOR_STEP;
                  decSpeed();
                  scheduleRedraw();
               }
               return TRUE;
            case REVERSED_ANIMATION:
               incSpeed();
               maxSpeed -= SPEED_INDICATOR_STEP;
               if (maxSpeed < -0.8)
                  maxSpeed = -MAX_SPEED_INDICATOR;
               scheduleRedraw();
               return TRUE;
            case PAUSED_ANIMATION:
               maxSpeed -= SPEED_INDICATOR_STEP;
               if (maxSpeed < -0.8)
                  maxSpeed = -MAX_SPEED_INDICATOR;
               if (beforePausing == REVERSED_ANIMATION) {
                  incSpeed();
               } else {
                  decSpeed();
                  if (animateBtwPtsPeriod >= MIN_SPEED)
                     beforePausing = REVERSED_ANIMATION;
               }
               scheduleRedraw();
               return TRUE;
            default:
               //fall through
               break;
            }
            break;

         case SoKeyboardEvent::E:
            this->escapeCallback(this->examinerObject);
            break;

         default:
            break; // To get rid of compiler warnings
         }
      }
      if (SoKeyboardEvent::isKeyReleaseEvent(ev, ke->getKey())) {
         switch (ke->getKey()) {
         case SoKeyboardEvent::LEFT_SHIFT:
            this->lshiftdown = false;
            return TRUE;
         case SoKeyboardEvent::RIGHT_SHIFT:
            this->rshiftdown = false;
            return TRUE;
         case SoKeyboardEvent::LEFT_CONTROL:
            this->lctrldown = false;
            return TRUE;
         case SoKeyboardEvent::RIGHT_CONTROL:
            this->rctrldown = false;
            return TRUE;
         default:
            break;
         }
      }
   }

   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == ROTATING)
      return FALSE;
   else
      return SoXtExaminerViewer::processSoEvent(ev);
}

// Called by hitting PageUp during animation.
void G4OpenInventorXtExaminerViewer::incSpeed() {
	if (std::ceil(animateBtwPtsPeriod * 100) >= 4) {
		if (speedStep > 0.08)
			speedStep -= 0.02;
		else
			speedStep = 0.02;
		animateBtwPtsPeriod -= speedStep;
	} else
		animateBtwPtsPeriod = 0.0;

	if (currentState != PAUSED_ANIMATION) {
		int lastIdx = refParticleTrajectory.size() - 1;
		if (refParticleIdx < lastIdx && !animateSensor->isScheduled())
			animateRefParticle();
	}
}

// Called by hitting PageDown during animation.
void G4OpenInventorXtExaminerViewer::decSpeed() {
	animateBtwPtsPeriod += speedStep;
	if (animateBtwPtsPeriod < MIN_SPEED) {
		if (std::floor(animateBtwPtsPeriod * 100) == 12) { // Errors in double representation
			speedStep = 0.08;
		} else if (animateBtwPtsPeriod > 0.12)
			speedStep += 0.02;
	} else {
		animateBtwPtsPeriod = MIN_SPEED;
		speedStep = START_STEP;
		maxSpeed = 0.0f;
		if (animateSensor->isScheduled())
			animateSensor->unschedule();
	}
}

// Based on the user's interaction the speed indicator bar needs to be adjusted.
void G4OpenInventorXtExaminerViewer::updateSpeedIndicator(void) {
	assert(this->sgeometry != NULL);

	SbVec3f * points = this->sgeometry->point.startEditing();

	if (points[10][0] == 0.0f)
		this->animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_ALL);
	if (points[14][0] == 0.0f)
		this->animSpeedSwitch->whichChild.setValue(SO_SWITCH_ALL);
	points[10][0] = this->maxSpeed;
	points[11][0] = this->maxSpeed;
	points[14][0] = this->maxSpeed;
	points[15][0] = this->maxSpeed;
	this->sgeometry->point.finishEditing();

	if (this->maxSpeed == 0.0f) {
		this->animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
		this->animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);
	}
}

void G4OpenInventorXtExaminerViewer::actualRedraw(void) {
	switch (currentState) {
	case ANIMATION:
	case REVERSED_ANIMATION:
	case PAUSED_ANIMATION:
		updateSpeedIndicator();
		SoXtExaminerViewer::actualRedraw();
		break;
	default:
		SoXtExaminerViewer::actualRedraw();
		break;
	}
}

void G4OpenInventorXtExaminerViewer::setReferencePath(SoLineSet *lineset, SoCoordinate3 *coords, bool append)
{
   // TODO:  Color the reference path
   // Disable the color stuff for now: changes all trajectories

// // We change the color of the trajectory too, so we get its material
//	nodeIndex = grpNode->findChild(trajectory);
//	SoMaterial * mat;
//	for(int i = 0; i < 100; ++i){
//		--nodeIndex;
//
//		tmpNode = grpNode->getChild(nodeIndex);
//		if(tmpNode->getTypeId() == SoMaterial::getClassTypeId()){
//			//node found
//			mat = (SoMaterial *)tmpNode;
//
//			break;
//		}
//	}
//
//
// // Restore default color for previously picked trajectory
// // and set different color for current pick
//	if(This->prevColorField)
//		((SoMFColor *)This->prevColorField)->setValue(0.0, 1.0, 0.0);
//	This->prevColorField = (void *)&mat->diffuseColor;
//
//
//	if(mat->diffuseColor.isConnected())
//		std::cout << "connected" << std::endl;
//
//	mat->diffuseColor.setValue(41.0/255.0, 230.0/255.0, 230.0/255.0);
//
//	std::cout << "R: " << mat->diffuseColor[0][0] << " ";
//	std::cout << "G: " << mat->diffuseColor[0][1] << " ";
//	std::cout << "B: " << mat->diffuseColor[0][2] << std::endl;

   // The trajectory is composed of all the polyline segments in the
   // multiple value field (SoMFInt32) numVertices.
   // For each of the numVertices.getNum()* polyline segments,
   // retrieve the points from the SoCoordinate3 node
   SbVec3f refParticlePt;

   if(!append)
      this->refParticleTrajectory.clear();

   for(int i = 0; i < lineset->numVertices.getNum(); ++i){
      for(int j = 0; j < lineset->numVertices[i]; ++j){
         refParticlePt = coords->point[j];
         this->refParticleTrajectory.push_back(refParticlePt);
      }
   }
   // Remove points that are too close to each other
   this->evenOutRefParticlePts();
   this->setReferencePathZPos();
   this->sortElements();
}


void G4OpenInventorXtExaminerViewer::setReferencePathZPos()
{
   refZPositions.clear();
   refZPositions.push_back(0);
   float dist;
   for(unsigned int i=0; i < this->refParticleTrajectory.size() - 1; ++i){
      dist = (refParticleTrajectory[i] - 
              refParticleTrajectory[i + 1]).length();
      refZPositions.push_back(refZPositions[i] + dist);
   }
}


void G4OpenInventorXtExaminerViewer::findAndSetRefPath()
{
   SoSearchAction action;
   action.setType(SoLineSet::getClassTypeId(),false);
   action.setInterest(SoSearchAction::ALL);
   action.apply(this->getSceneGraph());

   SoPathList &pathList = action.getPaths();

   if(pathList.getLength() != 0){

      SoCoordinate3 * coords = NULL;
      std::vector<SoCoordinate3 *> coordvec;
      std::vector<SoLineSet *> linevec;

      bool refPathFound = false;
      for(int i = 0; i < pathList.getLength(); ++i) {
         SoFullPath *path = (SoFullPath *)pathList[i];

         G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(path->getTail());
         for (size_t j = 0; j < attHolder->GetAttDefs().size(); ++j) {
            std::ostringstream oss;
            oss << G4AttCheck(attHolder->GetAttValues()[j],	attHolder->GetAttDefs()[j]);

            std::string findStr = "Type of trajectory (Type): ";
            std::string compareValue = "REFERENCE";
            size_t idx = oss.str().find(findStr);

            if(idx != std::string::npos) {
               if(oss.str().substr(idx + findStr.size(), compareValue.size()) == compareValue) {
                  coords = this->getCoordsNode(path);
                  if(coords != NULL){
                     refPathFound = true;
                     coordvec.push_back(coords);
                     linevec.push_back((SoLineSet *)path->getTail());
                  }
                  break;
               }
            }

            findStr = "Track ID (ID): ";
            idx = oss.str().find(findStr);
            if(idx != std::string::npos) {
               //index all primary tracks
               std::string tmpstr = oss.str().substr(idx + findStr.size(),1);
               std::istringstream buffer(tmpstr);
               int num;
               buffer >> num;
               if(num == 1) {

                  // Check if next character is a number, 
                  // in which case we don't have Track ID 1
                  // FWJ attempt to fix Coverity issue.
                  char nextChar = oss.str().at(idx+findStr.size()+1);
                  // const char * nextChar = 
                  // oss.str().substr(idx + findStr.size() + 1,1).c_str();
                  if(std::isdigit(nextChar))
                  // if(std::isdigit(nextChar[0]))
                     break;	//Not a primary track, continue with next track

                  coords = this->getCoordsNode(path);
                  if(coords != NULL){
                     coordvec.push_back(coords);
                     linevec.push_back((SoLineSet *)path->getTail());
                     break; //Found coords node, continue with next track
                  }
               }
               else
                  break;	//Not a primary track, continue with next track
            }
            else{
               //Not a Track ID attribute, fall through
            }
         }

         if(refPathFound)
            break;
      }

      if(coordvec.empty())
         return;		//No track with a Coordinate3 node found

      if(refPathFound){
         //set ref path to last traj, coord in the vecs
         this->setReferencePath(linevec.back(), coordvec.back());
         return;
      }
      //else

      int longestIdx = 0;
      float longestLength = 0.0;
      // For all paths
      for(unsigned int i=0;i < linevec.size(); ++i){

         //First generate a vector with all the points in this lineset
         std::vector<SbVec3f> trajectory;
         // For all lines in the i path
         for(int j=0; j < linevec[i]->numVertices.getNum(); ++j){
            // For all points in line j
            for(int k=0; k < linevec[i]->numVertices[j]; ++k){
               trajectory.push_back(coordvec[i]->point[k]);
            }
         }

         // Then calculate the total length
         float tmpLength=0.0;
         for(unsigned int j=0; j < trajectory.size() - 1; ++j){
            tmpLength += (trajectory[j] - trajectory[j + 1]).length();
         }

         if(tmpLength > longestLength){
            longestIdx = i;
            longestLength = tmpLength;
         }
      }

      // Set the longest path as the reference path
      this->setReferencePath(linevec[longestIdx], coordvec[longestIdx]);
   }
}


SoCoordinate3 * G4OpenInventorXtExaminerViewer::getCoordsNode(SoFullPath *path)
{
   SoLineSet *trajectory = (SoLineSet *)path->getTail();
   SoSeparator * grpNode = (SoSeparator*)(((SoFullPath*)path)->getNodeFromTail(1));
   int nodeIndex = grpNode->findChild(trajectory);
   SoNode * tmpNode;

   // We allow only 100 iterations, in case the node isn't found
   // (should take only a few iterations)
   for(int i = 0; i < 100; ++i){
      --nodeIndex;

      tmpNode = grpNode->getChild(nodeIndex);
      if(tmpNode->getTypeId() == SoCoordinate3::getClassTypeId()){
         //node found
         return (SoCoordinate3 *)tmpNode;
      }
   }
   return NULL;	//coords node not found
}


// Displays scene elements on the right side of listsDialog.
// else: scene graph is searched for Geant4_SoPolyhedron type nodes
void G4OpenInventorXtExaminerViewer::getSceneElements()
{
   std::string field, eltName;

   std::map<std::string, int> duplicates;
   std::map<std::string, int> sceneElts;
   SoSearchAction search;
   Geant4_SoPolyhedron *node;
   SoGroup *root = (SoGroup *)getSceneManager()->getSceneGraph();

   SoBaseKit::setSearchingChildren(TRUE);

   search.reset();
   search.setSearchingAll(TRUE);
   search.setInterest(SoSearchAction::ALL);
   search.setType(Geant4_SoPolyhedron::getClassTypeId(), 0);
   search.apply(root);

   SoPathList &pl = search.getPaths();


   // First find which names occur more than once so we can append a counter to them
   for(int i = 0; i < pl.getLength(); i++) {
      SoFullPath *path = (SoFullPath *)pl[i];
      node = (Geant4_SoPolyhedron *)path->getTail();
      eltName = node->getName();
      if(duplicates.count(eltName))
         duplicates[eltName]++;
      else
         duplicates[eltName] = 1;
   }

   for(int i = 0; i < pl.getLength(); i++) {
      float x,y,z;
      std::stringstream ssCount;
      SoFullPath *path = (SoFullPath *)pl[i];
      node = (Geant4_SoPolyhedron *)path->getTail();
      eltName = node->getName();
      field = eltName;
      if(duplicates[eltName] == 1)
         ssCount << "";//duplicates[field]
      else {
         if(sceneElts.count(eltName))
            sceneElts[eltName]++;
         else
            sceneElts[eltName] = 1;

         ssCount << sceneElts[eltName];
         field += "_";
      }

      field += ssCount.str();

      SoGetBoundingBoxAction bAction(getViewportRegion());
      bAction.apply(path);
      SbBox3f bBox = bAction.getBoundingBox();

      SbVec3f centr = bBox.getCenter();
      centr.getValue(x,y,z);

      path->ref();
      sceneElement el = { field, path, centr, 0.0 };
      this->sceneElements.push_back(el);
   }
}


float G4OpenInventorXtExaminerViewer::sqrlen(const SbVec3f &a)
{
   float x,y,z;
   a.getValue(x,y,z);
   return x*x + y*y + z*z;
}


void G4OpenInventorXtExaminerViewer::distanceToTrajectory(const SbVec3f &q,
                                                          float &dist,
                                                SbVec3f &closestPoint,
                                                          int &index)
{
   // a : Previous point on trajectory
   // b : Next point on trajectory
   // q : the point in space
   // dab, daq, dbq: distance between a & b, a & q, b & q
   //    
   // Theory:  A point p on a line ab is defined as:
   //
   // 				p(t) = a+t⋅(b–a)
   //
   // 			note: All are vectors except the parameter t
   //
   // When t is between 0 and 1 the point p is situated between a and b on ab.
   // The point p is defined in terms of the parameter t, subsequently so does
   // the distance from the query point q to the point p. To find the minimum
   // of that distance we differentiate it and set equal to zero:
   //
   //  			diff(Norm(p(t)- q)) = 0
   //
   //  		note: diff means taking the derivative with regard to t
   //
   // The resulting t is given in the code below. The square of the distance
   // between p and q is given by:
   //
   //  			d^2 = (Norm(p(t)-q))^2
   //
   // The expression found is given in the code below (current_dist)
   //
   // Ref: http://programmizm.sourceforge.net/blog/2012/
   //           distance-from-a-point-to-a-polyline
   //
   //    --PLG

   const size_t count = this->refParticleTrajectory.size();
   assert(count>0);

   SbVec3f b = this->refParticleTrajectory[0];
   SbVec3f dbq = b - q;
   float sqrDist = sqrlen(dbq);
   closestPoint = b;
   index = 0;
   for (size_t i = 1; i < count; ++i) {
      const SbVec3f a = b;
      const SbVec3f daq = dbq;
      b = this->refParticleTrajectory[i];
      dbq = b - q;
      const SbVec3f dab = a - b;

      float dab_x, dab_y, dab_z;
      dab.getValue(dab_x,dab_y,dab_z);
      float daq_x, daq_y, daq_z;
      daq.getValue(daq_x, daq_y, daq_z);
      float dbq_x, dbq_y, dbq_z;
      dbq.getValue(dbq_x, dbq_y, dbq_z);

      const float inv_sqrlen = 1./sqrlen(dab);
      const float t = (dab_x*daq_x + dab_y*daq_y + dab_z*daq_z)*inv_sqrlen;

      if (t<0.){
         // The trajectory point occurs before point a
         // Go to the next point
         continue;
      }
      float current_dist;
      if (t<=1.){
         // The trajectory point occurs between a and b.
         // Compute the distance to that point
         current_dist = daq_x*daq_x + daq_y*daq_y + daq_z*daq_z
            - t*(daq_x*dab_x + daq_y*dab_y + daq_z*dab_z)
            + t*t*(dab_x*dab_x + dab_y*dab_y + dab_z*dab_z);
      }
      else { //t>1.
         // The trajectory point occurs after b.
         // Get the distance to point b
         current_dist = sqrlen(dbq);
      }

      if (current_dist < sqrDist){
         sqrDist = current_dist;
         closestPoint = a + t*(b-a);
         index = i;
      }
   }

   dist = std::sqrt(sqrDist);
}


void G4OpenInventorXtExaminerViewer::sortElements()
{
   if(this->refParticleTrajectory.empty())
      return;

   float * trajLength = new float[this->refParticleTrajectory.size()];
   typedef std::map<elementForSorting, sceneElement> sortedMap;
   sortedMap sorted;

   // For every point on the reference trajectory, compute
   // the total length from the start
   SbVec3f prevPoint;
   std::vector<SbVec3f>::iterator itRef = this->refParticleTrajectory.begin();
   int trajIndex = 0;
   prevPoint = *itRef;
   trajLength[trajIndex] = 0.0;
   ++itRef;
   ++trajIndex;
   for(; itRef != this->refParticleTrajectory.end(); ++itRef, ++trajIndex){
      trajLength[trajIndex] = trajLength[trajIndex-1] + (*itRef - prevPoint).length();
      prevPoint = *itRef;
   }

   // Compute the smallest distance between the element
   // and the reference trajectory (find the closest point),
   // then map the element to the trajectory length of that
   // point (calculated above)
   SoGetBoundingBoxAction bAction(this->getViewportRegion());
   SbVec3f elementCoord;
   std::vector<sceneElement>::iterator itEl;
   int elementIndex;
   elementForSorting el;
   for(itEl = this->sceneElements.begin(), elementIndex = 0;
       itEl != this->sceneElements.end(); ++itEl, ++elementIndex){
      bAction.apply(itEl->path);

      // FWJ sceneElement already has a center
      elementCoord = itEl->center;
      // ... and this sometimes returns an empty box!
      //      elementCoord = bAction.getBoundingBox().getCenter();
      //      if (bAction.getBoundingBox().isEmpty()) {
      //         G4cout << "sortElements: Box is empty!" << G4endl;
      //         G4cout << "   element name=" << itEl->name << G4endl;
      //      }

      int index;
      distanceToTrajectory(elementCoord, el.smallestDistance, el.closestPoint, index);
      itEl->closestPointZCoord = el.closestPointZCoord = trajLength[index];
      el.distanceToBeamlineStart = (itEl->center - this->refParticleTrajectory[0]).length();

      // This map of the scene elements (or their coordinates rather)
      // is automatically sorted by trajectory length (Z coord), then
      // by the distance between the element and the point in case the Z coord
      // is the same as another element.  This is done by using as a key
      // an element structure which implements the operator for weak ordering
      sorted.insert(std::make_pair(el,*itEl));
   }

   // store the sorted elements into the vector field
   this->sceneElements.clear();

   sortedMap::iterator itSorted = sorted.begin();
   for(; itSorted != sorted.end(); itSorted++)
      this->sceneElements.push_back(itSorted->second);

   this->zcoordSetFlag = true;


   Widget formTop = XtNameToWidget(this->listsDialog, "FormTop");
   Widget formTopRight = XtNameToWidget(formTop, "FormTopRight");

   this->createElementsList(formTopRight);

   delete[] trajLength;
}


void G4OpenInventorXtExaminerViewer::createElementsList(Widget formTopRight)
{
   if(this->myElementList != NULL)
      XtUnmanageChild(this->myElementList);

   int size = this->sceneElements.size();
   XmString *elements = (XmString *) XtMalloc(size * sizeof(XmString));

   std::vector<sceneElement>::const_iterator it;
   int count = 0;
   std::stringstream ss;
   for(it=this->sceneElements.begin(); it!=this->sceneElements.end(); ++it) {
      ss << it->name;
      if(zcoordSetFlag)
         ss << " [" << it->closestPointZCoord << "]";
      elements[count] = XmStringCreateLocalized((char *)ss.str().c_str());
      ++count;
      ss.str("");
   }

   Arg args[10];
   int n;

   // Label Right
   n = 0;
   Widget labelRight;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM);	n++;

   labelRight = XmCreateLabelGadget(formTopRight, (char*)"Element [S mm]",
                                    args, n);
   XtManageChild(labelRight);

   // List Right
   n = 0;
   XtSetArg(args[n], XmNvisibleItemCount, 7);	n++;
   XtSetArg(args[n], XmNitemCount, size);	n++;
   XtSetArg(args[n], XmNitems, elements);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNtopWidget, labelRight);	n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM);	n++;
   // FWJ
   XtSetArg(args[n], XmNwidth, 240);	n++;
   //   XtSetArg(args[n], XmNwidth, 280);	n++;
   //	XtSetArg(args[n], XmNwidth, 300);	n++;

   this->myElementList = XmCreateScrolledList(formTopRight, (char *) "ListRight", args, n);

   XtAddCallback(this->myElementList, XmNbrowseSelectionCallback,
                 (XtCallbackProc) lookAtSceneElementCB, this);
   xmAddMouseEventHandler(this->myElementList); // Add scrolling functionality
   XtManageChild(this->myElementList);

   if (elements != NULL) {
      for (int i = 0; i < size; i++)
         XmStringFree(elements[i]);
      XtFree((char *) elements);
   }
}


// Pops up a custom dialog listsDialog containing 
// scene elements and viewpoints.

void G4OpenInventorXtExaminerViewer::constructListsDialog(Widget w,
                                             XtPointer client_data,
                                             XtPointer)
{
   // G4cout << "DEBUG constructListsDialog w = " << w << G4endl;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   if (This->listsDialog) {
      return;
   }

   if (This->currentState == ANIMATION || This->currentState == PAUSED_ANIMATION) {
      if (This->animateSensor->isScheduled())
         This->animateSensor->unschedule();
      This->refParticleIdx = This->prevRefIdx;
      This->restoreCamera();
      This->currentState = This->prevState;
   }

   This->step = 1; // Default values
   This->refParticleIdx = 0;
   if (This->refParticleTrajectory.size()){
      This->prevPt = This->refParticleTrajectory[0]; // For calculating distance
   }

   This->getSceneElements();

   int n = 0;
   Arg args[10];
   Atom WM_DELETE_WINDOW;

   ///////////////////////CUSTOM listsDialog///////////////////////////////

   Widget topShell;
   // FWJ gets the topmost window containing This->getParentWidget()
   // This is unnecessary because the parent is passed in
   //   topShell = SoXt::getShellWidget(This->getParentWidget());
   topShell = w;
   // G4cout << "DEBUG PARENT (topShell) FOR AUX WINDOW = " << topShell << G4endl;

   // Shell Dialog
   std::string dialogNameStr = This->fileName.substr(This->fileName.rfind('/') + 1);
   const int nDialog = dialogNameStr.size() + 1;
   char *dialogName = new char[nDialog];
   strncpy(dialogName, dialogNameStr.c_str(), nDialog);

   n = 0;
   XtSetArg(args[n], XmNx, 610);	n++;
   This->myShellDialog = XmCreateDialogShell(topShell, dialogName, args, n);

   delete[] dialogName;
   WM_DELETE_WINDOW = XInternAtom(XtDisplay(w), "WM_DELETE_WINDOW", False);
   XmAddWMProtocolCallback(This->myShellDialog, WM_DELETE_WINDOW,
                           (XtCallbackProc)closeListsDialogCB, This);

   // Main Pane(listsDialog)
   n = 0;
   XtSetArg(args[n], XmNsashWidth, 1);	n++;
   XtSetArg(args[n], XmNsashHeight, 1);	n++;
   XtSetArg(args[n], XmNseparatorOn, False);	n++;
   // FWJ
   This->listsDialog = XmCreatePanedWindow(This->myShellDialog, (char *) "MainPane",
                                           args, n);


   ////////////////////////TOP FORM//////////////////////////
   n = 0;
   // FWJ fails compile
   //   Widget formTop = XmCreateForm(This, (char *) "FormTop", args, n);
   Widget formTop = XmCreateForm(This->listsDialog, (char *) "FormTop", args, n);

   n = 0;
   XtSetArg(args[n], XmNmarginWidth, 8);	n++;
   XtSetArg(args[n], XmNmarginHeight, 8);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM);	n++;
   Widget formTopRight = XmCreateForm(formTop, (char *) "FormTopRight", args,
                                      n);

   n = 0;
   XtSetArg(args[n], XmNmarginWidth, 8);	n++;
   XtSetArg(args[n], XmNmarginHeight, 8);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNrightWidget, formTopRight);	n++;
   XtSetArg(args[n], XmNrightOffset, 10);	n++;
   XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM);	n++;
   Widget formTopLeft = XmCreateForm(formTop, (char *) "FormTopLeft", args, n);

   /////TOP RIGHT/////

   This->createElementsList(formTopRight);
   XtManageChild(formTopRight);

   /////TOP LEFT/////

   // Label Left
   n = 0;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM);	n++;
   Widget labelLeft = XmCreateLabelGadget(formTopLeft, (char *) "ViewPoints",
                                          args, n);
   XtManageChild(labelLeft);

   // List Left
   n = 0;
   XtSetArg(args[n], XmNlistSizePolicy, XmRESIZE_IF_POSSIBLE);	n++;
   XtSetArg(args[n], XmNvisibleItemCount, 7);	n++;
   //	XtSetArg(args[n], XmNwidth, 140);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNtopWidget, labelLeft);	n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNrightWidget, This->myElementList);	n++;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM);	n++;
   // FWJ
   XtSetArg(args[n], XmNwidth, 160);	n++;
   // XtSetArg(args[n], XmNwidth, 200);	n++;

   This->myViewPtList = XmCreateScrolledList(formTopLeft, (char *) "ListLeft",
                                             args, n);
   if (This->viewPtList.size())
      This->addViewPoints();
   XtAddCallback(This->myViewPtList, XmNbrowseSelectionCallback,
                 (XtCallbackProc) loadBookmarkCB, This);
   xmAddMouseEventHandler(This->myViewPtList); // Add scrolling functionality

   XtManageChild(This->myViewPtList);

   XtManageChild(formTopLeft);

   XtManageChild(formTop);

   ////////////////////MIDDLE FORM///////////////////////////
   n = 0;
   XtSetArg(args[n], XmNmarginWidth, 6);	n++;
   // FWJ fails compile
   //   Widget formMiddle = XmCreateForm(This->canvas, (char *) "MiddleForm", args, n);
   Widget formMiddle = XmCreateForm(This->listsDialog, (char *) "MiddleForm", args, n);

   // Label
   n = 0;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNtopWidget, This->myViewPtList);	n++;
   Widget label = XmCreateLabelGadget(formMiddle, (char *) "Selection", args,
                                      n);
   XtManageChild(label);

   // Text
   n = 0;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNtopWidget, label);	n++;
   XtSetArg(args[n], XmNtopOffset, 3);	n++;
   XtSetArg(args[n], XmNmaxLength, This->MAX_VP_NAME);	n++;
   This->viewPtSelection = XmCreateText(formMiddle, (char *) "Txt", args, n);
   XtManageChild(This->viewPtSelection);

   Dimension h1, h2, h;
   XtVaGetValues(label, XmNheight, &h1, NULL);
   XtVaGetValues(This->viewPtSelection, XmNheight, &h2, NULL);

   h = (Dimension) (1.1 * (h1 + h2));

   XtVaSetValues(formMiddle, XmNpaneMaximum, h, XmNpaneMinimum, h, NULL);
   XtManageChild(formMiddle);

   /////////////////////BOTTOM FORM///////////////////////////
   // Action Area Form
   n = 0;
   XtSetArg(args[n], XmNfractionBase, 4);	n++;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_WIDGET);	n++;
   XtSetArg(args[n], XmNtopWidget, This->viewPtSelection);	n++;
   // FWJ fails compile
   //   Widget formAction = XmCreateForm(This, (char *) "ActionForm", args, n);
   Widget formAction = XmCreateForm(This->listsDialog, (char *) "ActionForm", args, n);

   n = 0;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM);	n++;
   XtSetArg(args[n], XmNtopOffset, 3);	n++;
   XtSetArg(args[n], XmNbottomOffset, 5);	n++;
   Widget separator = XmCreateSeparatorGadget(formAction, (char *) "Sep", args, n);

   XtManageChild(separator);

   Widget button = XmCreatePushButton(formAction, (char *) "Delete", NULL, 0);
   XtVaSetValues(button, XmNtopAttachment, XmATTACH_WIDGET, XmNtopWidget,
                 separator, XmNbottomAttachment, XmATTACH_FORM, XmNleftAttachment,
                 XmATTACH_POSITION, XmNleftPosition, 0, XmNrightAttachment,
                 XmATTACH_POSITION, XmNrightPosition, 1,
                 XmNdefaultButtonShadowThickness, 2, XmNwidth, 40, XmNheight, 30,
                 NULL);

   XtAddCallback(button, XmNactivateCallback,
                 (XtCallbackProc) deleteBookmarkCB, This);
   XtManageChild(button);

   button = XmCreatePushButton(formAction, (char *) "Rename", NULL, 0);
   XtVaSetValues(button, XmNtopAttachment, XmATTACH_WIDGET, XmNtopWidget,
                 separator, XmNbottomAttachment, XmATTACH_FORM, XmNleftAttachment,
                 XmATTACH_POSITION, XmNleftPosition, 1, XmNrightAttachment,
                 XmATTACH_POSITION, XmNrightPosition, 2,
                 XmNdefaultButtonShadowThickness, 2, XmNwidth, 40, XmNheight, 30,
                 NULL);

   XtAddCallback(button, XmNactivateCallback,
                 (XtCallbackProc) renameBookmarkCB, This);
   XtManageChild(button);

   button = XmCreatePushButton(formAction, (char *) "Sort", NULL, 0);
   XtVaSetValues(button, XmNtopAttachment, XmATTACH_WIDGET, XmNtopWidget,
                 separator, XmNbottomAttachment, XmATTACH_FORM, XmNleftAttachment,
                 XmATTACH_POSITION, XmNleftPosition, 2, XmNrightAttachment,
                 XmATTACH_POSITION, XmNrightPosition, 3,
                 XmNdefaultButtonShadowThickness, 2, XmNwidth, 40, XmNheight, 30,
                 NULL);

   XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) sortBookmarksCB, This);
   XtManageChild(button);

   button = XmCreatePushButton(formAction, (char *) "Close", NULL, 0);
   XtVaSetValues(button, XmNtopAttachment, XmATTACH_WIDGET, XmNtopWidget,
                 separator, XmNbottomAttachment, XmATTACH_FORM, XmNleftAttachment,
                 XmATTACH_POSITION, XmNleftPosition, 3, XmNrightAttachment,
                 XmATTACH_POSITION, XmNrightPosition, 4,
                 XmNdefaultButtonShadowThickness, 2, XmNwidth, 40, XmNheight, 30,
                 NULL);

   XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) closeListsDialogCB, This);
   XtManageChild(button);

   XtManageChild(formAction);
   XtVaGetValues(button, XmNheight, &h1, NULL);
   XtVaSetValues(formAction, XmNpaneMaximum, h1, XmNpaneMinimum, h1, NULL);

   XtManageChild(This->listsDialog);

   ////////////////////////CUSTOM listsDialog///////////////////////////////
}


// Called when user clicks a scene element in listsDialog.
// Zooms onto that element.
void G4OpenInventorXtExaminerViewer::lookAtSceneElementCB(Widget,
                                           XtPointer client_data,
                                           XtPointer call_data)
{
   char *value;
   std::string elementField;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   SoCamera * cam = This->getCamera();

   if (This->SoXtExaminerViewer::isAnimating())
      This->stopAnimating();

   XmListCallbackStruct *cbs = (XmListCallbackStruct *) call_data;

   value = (char *) XmStringUnparse(cbs->item, XmFONTLIST_DEFAULT_TAG,
                                    XmCHARSET_TEXT, XmCHARSET_TEXT, NULL, 0, XmOUTPUT_ALL);
   if (This->currentState == ANIMATION || This->currentState == REVERSED_ANIMATION
       || This->currentState == PAUSED_ANIMATION ) {
      if (This->animateSensor->isScheduled())
         This->animateSensor->unschedule();
      This->setSuperimpositionEnabled(This->superimposition, FALSE);
      This->maxSpeed = 0.0f;
      This->scheduleRedraw();
      This->restoreCamera();
      This->currentState = This->prevState;
   } else if (This->currentState == VIEWPOINT)
      This->setSuperimpositionEnabled(This->superimposition, FALSE);

   elementField = value;

   int idx = elementField.find_last_of("[");
   if(idx == -1)
      idx = elementField.size(); //if "[" not found for whatever reason (list not sorted)
   else
      idx--; // To get rid of the space that is between the name and '['

   bool error = false;
   SoFullPath *path;
   SoSearchAction search;
   SoNode *root = This->getSceneManager()->getSceneGraph();
   int counter, idxUnderscore = elementField.find_last_of("_");

   This->parseString<int>(counter, elementField.substr(idxUnderscore + 1, idx), error);

   SoBaseKit::setSearchingChildren(TRUE);
   search.reset();
   search.setSearchingAll(TRUE);

   if(error) { // No counter is present => element name was not modified
      This->curEltName = elementField.substr(0, idx);
      search.setName(This->curEltName.c_str());
      search.apply(root);

      path = (SoFullPath *)search.getPath();
   }
   else {
      This->curEltName = elementField.substr(0, idxUnderscore);
      search.setInterest(SoSearchAction::ALL);
      search.setName(This->curEltName.c_str());
      search.apply(root);

      SoPathList &pl = search.getPaths();
      path = (SoFullPath *)pl[counter - 1]; // Since counter starts at 1, not 0
   }

   G4ThreeVector global;

   if ((idx > 0) && (path)) {

      if(!This->refParticleTrajectory.empty()){

         SoGetBoundingBoxAction bAction(This->getViewportRegion());
         bAction.apply(path);
         SbBox3f bBox = bAction.getBoundingBox();
         SbVec3f elementCoord = bBox.getCenter();

         This->refParticleIdx = 0;
         SbVec3f p;

         float absLengthNow, absLengthMin;
         int maxIdx = This->refParticleTrajectory.size() - 2;
         int targetIdx = 0;
         SbVec3f dir;

         p = This->refParticleTrajectory[This->refParticleIdx];
         absLengthMin = (p - elementCoord).length();
         This->refParticleIdx++;

         // Find a ref. particle's point closest to element's global coords
         while (This->refParticleIdx < maxIdx) {
            p = This->refParticleTrajectory[This->refParticleIdx];
            absLengthNow = (p - elementCoord).length();

            if (absLengthNow < absLengthMin) {
               absLengthMin = absLengthNow;
               targetIdx = This->refParticleIdx;
            }
            This->refParticleIdx++;
         }

         if (This->currentState != BEAMLINE) { // Set up default zoom
            SbVec3f p1, pN;
            This->currentState = BEAMLINE;
            This->prevParticleDir = SbVec3f(0,0,0); //so that moveCamera() knows sets default parameters
            
            p1 = This->prevPt = This->refParticleTrajectory[0];
            pN = This->refParticleTrajectory[This->refParticleTrajectory.size() - 1];
            This->distance = (pN - p1).length() / 10;

            // FWJ Rather than switching to a default height, it is more flexible
            // to keep the same height(magnification) while moving the camera.
            // if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
            //    ((SoOrthographicCamera *) cam)->height.setValue(This->defaultHeight);
            // // FWJ Restore the default height instead of hard-wired value
            // // ((SoOrthographicCamera *) cam)->height.setValue(10000.0f);
            // }
            // else if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))

            // FWJ required to avoid extreme perspective after camera move:
            if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))
               ((SoPerspectiveCamera*)cam)->heightAngle.setValue(This->defaultHeightAngle);

         } else {
            if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))
               This->distance = (This->prevPt - cam->position.getValue()).length();
         }
         This->refParticleIdx = targetIdx;

         //////////////////////////////////////////////////////////////
         This->setSuperimpositionEnabled(This->superimposition, TRUE);
         This->axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
         This->animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
         This->animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);
         This->scheduleRedraw();
         //////////////////////////////////////////////////////////////

         This->moveCamera(This->distance);
         XtFree(value);

      }
      
      else{
         This->offsetFromCenter.setValue(0, 0, 1);
         This->distance = 50;// small number since using viewAll() for default zoom
         This->upVector.setValue(0, 1, 0);
         This->moveCamera(This->distance);
         cam->viewAll(path, This->getViewportRegion());
      }
   }

   XmTextSetString(This->viewPtSelection, NULL);
}


// Destroyes listsDialog and resets necessary member fields.

void G4OpenInventorXtExaminerViewer::closeListsDialogCB(Widget, 
                                         XtPointer client_data,
                                         XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   This->sceneElements.clear();
   This->refParticleTrajectory.clear();

   This->currentState = GENERAL;
   XtDestroyWidget(This->myShellDialog);
   This->listsDialog = NULL;
}

// Called when user clicks left arrow button. Loads previous viewpoint.
void G4OpenInventorXtExaminerViewer::prevViewPtCB(Widget, XtPointer client_data,
                                                  XtPointer) {
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (This->viewPtIdx == 0)
      This->viewPtIdx = This->viewPtList.size() - 1;
   else
      This->viewPtIdx--;

   This->writeViewPtIdx();
   This->setViewPt();
}

// Called when user clicks right arrow button. Loads next viewpoint.
void G4OpenInventorXtExaminerViewer::nextViewPtCB(Widget, XtPointer client_data,
                                                  XtPointer) {
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (This->viewPtIdx >= (int) This->viewPtList.size() - 1)
      This->viewPtIdx = 0;
   else
      This->viewPtIdx++;

   This->writeViewPtIdx();
   This->setViewPt();
}


// Updates the viewPtIdx in a viewpoint file.

void G4OpenInventorXtExaminerViewer::writeViewPtIdx()
{
   std::string idxStr;
   std::stringstream out;
   out << viewPtIdx;
   idxStr = out.str();
   fileOut.seekp(0, std::ios::beg);

   while ((int) idxStr.length() < MAX_VP_IDX) {
      idxStr += " ";
   }

   fileOut << idxStr << "\n";
   fileOut.flush();
   fileOut.seekp(0, std::ios::end);
}


// Sets the viewpoint based on camera data that viewPtIdx is pointing to.

void G4OpenInventorXtExaminerViewer::setViewPt()
{
   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == ROTATING) {

      if (animateSensor->isScheduled())
         animateSensor->unschedule();
      setSuperimpositionEnabled(superimposition, FALSE);
      maxSpeed = 0.0f;
      scheduleRedraw();
   }

   SoCamera * camera = getCamera();
   if (camera == NULL) {
      String dialogName = (char *) "Missing Camera Node";
      std::string msg = "Camera is null. Unable to set the viewpoint.";
      warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   if (!viewPtList.size()) {
      String dialogName = (char *) "Missing Viewpoints";
      std::string msg = "There are no viewpoints to load.";
      warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   if (SoXtExaminerViewer::isAnimating())
      stopAnimating();

   if (currentState != VIEWPOINT) {
      currentState = VIEWPOINT;
      //////////////////////////////////////////////////////////////
      setSuperimpositionEnabled(superimposition, TRUE);
      axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
      animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
      animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);

      scheduleRedraw();
      ///////////////////////////////////////////////////////////////
   }

   curViewPtName = viewPtList[viewPtIdx].viewPtName;
   camera->viewportMapping = viewPtList[viewPtIdx].viewportMapping;
   camera->position = viewPtList[viewPtIdx].position;
   camera->orientation = viewPtList[viewPtIdx].orientation;
   camera->aspectRatio = viewPtList[viewPtIdx].aspectRatio;
   camera->nearDistance = viewPtList[viewPtIdx].nearDistance;
   camera->farDistance = viewPtList[viewPtIdx].farDistance;
   camera->focalDistance = viewPtList[viewPtIdx].focalDistance;

   // Restore camera height (changed by zooming)
   if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
      if (viewPtList[viewPtIdx].camType == ORTHOGRAPHIC) {
         toggleCameraType();
         camera = getCamera();
         ((SoOrthographicCamera *) camera)->height.setValue(
                                                            viewPtList[viewPtIdx].height);
      } else
         ((SoPerspectiveCamera *) camera)->heightAngle.setValue(
                                                                viewPtList[viewPtIdx].height);
   } else if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
      if (viewPtList[viewPtIdx].camType == PERSPECTIVE) {
         toggleCameraType();
         camera = getCamera();
         ((SoPerspectiveCamera *) camera)->heightAngle.setValue(
                                                                viewPtList[viewPtIdx].height);
      } else
         ((SoOrthographicCamera *) camera)->height.setValue(
                                                            viewPtList[viewPtIdx].height);
   } else {
      SoDebugError::post("G4OpenInventorXtExaminerViewer::setViewPt",
                         "Only Perspective and Orthographic cameras are supported.");
      return;
   }

}


// Pops up a prompt asking for a new viewpoint name.

void G4OpenInventorXtExaminerViewer::saveViewPtCB(Widget w, 
                                     XtPointer client_data,
                                     XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (This->fileName.empty()) {
      newViewPtFileCB(w, This, NULL);
      This->returnToSaveVP = true;
      return; // Need to return and call this fn again from newViewPtFileCB since flow of control does not stall here but keeps going
   }

   int n = 0;
   Arg args[4];
   Widget nameViewPtDialog;
   Widget parent = This->getParentWidget();	//gets the dialogshell of the ExaminerViewer widget
   XmString label = XmStringCreateLocalized((char *) "Name the viewpoint:");

   XtSetArg(args[n], XmNselectionLabelString, label);	n++;
// Prevent the dialog from closing automatically, in case the name is wrong
   XtSetArg(args[n], XmNautoUnmanage, False);	n++;
// FWJ
   XtSetArg(args[n], XmNtitle, "Save Bookmark"); n++;
   nameViewPtDialog = XmCreatePromptDialog(parent, String("Save Bookmark"),
                                           args, n);

   XmStringFree(label);
   XtAddCallback(nameViewPtDialog, XmNokCallback, getViewPtNameCB, This);
   XtAddCallback(nameViewPtDialog, XmNcancelCallback,
                 (XtCallbackProc) XtDestroyWidget, NULL);

   Widget text = XtNameToWidget(nameViewPtDialog, "Text");
   XtVaSetValues(text, XmNmaxLength, This->MAX_VP_NAME, NULL);
   std::string autoName = "";
   if (!This->warningFlag) { //leave the TextField as it is if coming back from warning dialog
      autoName = This->viewPtAutoName();
   }
   This->warningFlag = false;
   XmTextSetString(text, (char *) autoName.c_str());
   XmTextSetInsertionPosition(text, autoName.length());

   XtUnmanageChild(XtNameToWidget(nameViewPtDialog, "Help"));
   XtManageChild(nameViewPtDialog);
}


std::string G4OpenInventorXtExaminerViewer::viewPtAutoName()
{
   std::string viewPt;
   std::stringstream sstream;
   std::vector<int> existingViewPts;
   int tmp;

   //Build the list of names of the form viewpoint_* already present
   for (unsigned int i = 0; i < this->viewPtList.size(); ++i) {
      viewPt = this->viewPtList[i].viewPtName;
      if (viewPt.find("viewpoint_") != std::string::npos) {
         tmp = atoi(viewPt.substr(10).c_str());
         if (tmp == 0) {
            //0 means couldn't convert to integer OR viewpoint_0
            if (!viewPt.compare("viewpoint_0"))
               existingViewPts.push_back(0);
         } else
            existingViewPts.push_back(tmp);
      }
   }

   sstream.str("");
   sstream.clear();

   //Return the view viewpoint_* name available
   if (existingViewPts.size() > 0) {
      int vpNum = 0;
      while (true) {
         if (std::find(existingViewPts.begin(), existingViewPts.end(), vpNum)
             == existingViewPts.end()) {
            sstream << "viewpoint_" << vpNum;
            return sstream.str();
         }
         ++vpNum;
      }
   } else {
      return "viewpoint_0";
   }
   return "";
}


void G4OpenInventorXtExaminerViewer::abbrOutputCB(Widget,
                                   XtPointer client_data,
                                   XtPointer)
{
   G4OpenInventorXtExaminerViewer * This =
      (G4OpenInventorXtExaminerViewer *) client_data;
// G4cout << "DISARMCALLBACK abbrOutputFlag=" << This->abbrOutputFlag << G4endl;
   This->abbrOutputFlag = !(This->abbrOutputFlag);
}


void G4OpenInventorXtExaminerViewer::pickRefPathCB(Widget,
                                    XtPointer client_data,
                                    XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   // Save viewing state and go to picking mode
   This->viewingBeforePickRef = This->isViewing();
   if(This->isViewing())
      This->setViewing(false);
   This->setComponentCursor(SoXtCursor(SoXtCursor::CROSSHAIR));
   This->pickRefPathFlag = true;
}


void G4OpenInventorXtExaminerViewer::switchWireFrameCB(Widget w,
                                    XtPointer client_data,
                                    XtPointer)
{
   G4OpenInventorXtExaminerViewer* This = 
      (G4OpenInventorXtExaminerViewer*)client_data;
   //   xmToggleButton theToggleButton = (xmToggleButton)w;
   if (XmToggleButtonGetState(w)) {
         This->setDrawStyle(SoXtViewer::STILL, SoXtViewer::VIEW_LINE);
         This->setDrawStyle(SoXtViewer::INTERACTIVE, SoXtViewer::VIEW_LINE);
      } else {
         This->setDrawStyle(SoXtViewer::STILL, SoXtViewer::VIEW_AS_IS);
         This->setDrawStyle(SoXtViewer::INTERACTIVE,
                            SoXtViewer::VIEW_SAME_AS_STILL);
      }
}


// Examines new viewpoint name and if OK calls saveViewPt.

void G4OpenInventorXtExaminerViewer::getViewPtNameCB(Widget w, 
                                        XtPointer client_data,
                                        XtPointer call_data)
{
   char *name = NULL;
   std::string strName;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   XmSelectionBoxCallbackStruct *cbs =
      (XmSelectionBoxCallbackStruct *) call_data;
   XmStringGetLtoR(cbs->value, XmFONTLIST_DEFAULT_TAG, &name);

   if (!name) {
      return;
   }
   if (!*name) {
      XtFree(name);
      return;
   }

   strName = name;
   XtFree(name);

   int beg = strName.find_first_not_of(' '); // Remove leading/trailing spaces
   int end = strName.find_last_not_of(' ');
   strName = strName.substr(beg, end - beg + 1);

   bool nameExists = false;
   int size = This->viewPtList.size();
   for (int i = 0; i < size; i++) {
      if (!strcmp(This->viewPtList[i].viewPtName, strName.c_str())) {
         nameExists = true;
         break;
      }
   }

   if (!nameExists) {
      const int nVPName = This->MAX_VP_NAME + 1;
      name = new char[nVPName];
      strncpy(name, strName.c_str(), nVPName);
      if (This->viewPtIdx == -1)
         This->viewPtIdx = 0;
      This->saveViewPt(name);
      if (This->listsDialog) {
         XmListAddItemUnselected(This->myViewPtList, cbs->value, 0); // vpName
      }
      //Dismiss the nameViewPtDialog dialog
      XtUnmanageChild(w);
   } else {
      String dialogName = (char *) "Existing Viewpoint";
      std::string msg = "The viewpoint already exists.";
      This->warningMsgDialog(msg, dialogName, NULL);

   }
}


// Saves current camera parameters to a viewpoint file.

void G4OpenInventorXtExaminerViewer::saveViewPt(char *name)
{
   SbVec3f axis;
   viewPtData tmp;
   float x, y, z, angle;
   SoCamera * camera = getCamera();

   if (viewPtList.size() == 0) {
      writeViewPtIdx();
      XtSetSensitive(nextViewPtButton, True); // Makes arrow buttons clickable
      XtSetSensitive(prevViewPtButton, True);
   }

   tmp.viewPtName = name;
   tmp.viewportMapping = camera->viewportMapping.getValue();
   tmp.position = camera->position.getValue();
   tmp.orientation = camera->orientation.getValue();
   tmp.aspectRatio = camera->aspectRatio.getValue();
   tmp.nearDistance = camera->nearDistance.getValue();
   tmp.farDistance = camera->farDistance.getValue();
   tmp.focalDistance = camera->focalDistance.getValue();

   // Save camera height (changed by zooming)
   if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
      tmp.height = ((SoPerspectiveCamera *) camera)->heightAngle.getValue();
      tmp.camType = PERSPECTIVE;
   } else if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
      tmp.height = ((SoOrthographicCamera *) camera)->height.getValue();
      tmp.camType = ORTHOGRAPHIC;
   } else {
      SoDebugError::post("G4OpenInventorXtExaminerViewer::saveViewPtCB",
                         "Only Perspective and Orthographic cameras are supported.");
      return;
   }

   viewPtList.push_back(tmp);

   // Now save the view point to a .txt file
   std::string vpName = name;

   while ((int) vpName.size() <= MAX_VP_NAME)
      vpName += " ";

   fileOut << vpName << std::endl;
   tmp.position.getValue(x, y, z);
   fileOut << x << " " << y << " " << z << std::endl;

   // Reusing x, y and z for storing the axis
   tmp.orientation.getValue(axis, angle);
   axis.getValue(x, y, z);
   fileOut << x << " " << y << " " << z << " " << angle << std::endl;

   fileOut << tmp.camType << " " << tmp.height << std::endl;
   fileOut << tmp.focalDistance << " ";
   fileOut << tmp.nearDistance << " ";
   fileOut << tmp.farDistance << std::endl;
   fileOut << tmp.viewportMapping << " ";
   fileOut << tmp.aspectRatio << "\n" << std::endl;
   fileOut.flush();
   viewPtIdx++;
}


void G4OpenInventorXtExaminerViewer::deleteViewPtCB(Widget,
                                     XtPointer client_data,
                                     XtPointer)
{
   G4OpenInventorXtExaminerViewer * This =
      (G4OpenInventorXtExaminerViewer *) client_data;
   This->deleteViewPt();
}


// Deletes current viewpoint the user is looking at.
// Updates the input file and bookmarks as well.

void G4OpenInventorXtExaminerViewer::deleteViewPt(char *vpName)
{
   std::string line;
   int end;
   fileIn.open(fileName.c_str());
   std::ofstream out("temporaryFile.txt");

   if (!vpName)
      vpName = viewPtList[viewPtIdx].viewPtName;

   if (listsDialog) {
      XmString vpNameStr = XmStringCreateLocalized(vpName);

      XmListDeleteItem(myViewPtList, vpNameStr);
      XmStringFree(vpNameStr);
   }

   getline(fileIn, line); // Printing the viewpoint idx
   out << line << "\n";

   while (getline(fileIn, line)) {
      end = line.find_last_not_of(' ');
      line = line.substr(0, end + 1);
      if (!strcmp(line.c_str(), vpName)) { // Equal
         while (line.size()) {
            getline(fileIn, line);
         }

         while (getline(fileIn, line))
            out << line << "\n";
      } else {
         while (line.size()) {
            out << line << "\n";
            getline(fileIn, line);
         }
         out << "\n";
      }
   }

   int idx = 0; // Remove viewpoint from the vector
   int size = viewPtList.size();
   while (idx < size) {
      if (!strcmp(viewPtList[idx].viewPtName, vpName)) {
         viewPtList.erase(viewPtList.begin() + idx);
         break;
      }
      idx++;
   }

   out.close();
   fileOut.close();
   fileIn.clear();
   fileIn.close();

   // FWJ check return status
   int istat = remove(fileName.c_str());
   if (istat == -1) {
      char dialogName[] = "Warning";
      warningMsgDialog("Error removing bookmarks file", dialogName,
                       NULL);
   }
   istat = rename("temporaryFile.txt", fileName.c_str());
   if (istat == -1) {
      char dialogName[] = "Warning";
      warningMsgDialog("Error renaming bookmarks file", dialogName,
                       NULL);
   }
   fileOut.open(fileName.c_str(), std::ios::in);
   fileOut.seekp(0, std::ios::end);

   if (!viewPtList.size()) { // viewPtList is empty
      curViewPtName = (char *) "";
      scheduleRedraw();
      XtSetSensitive(nextViewPtButton, False);
      XtSetSensitive(prevViewPtButton, False);
   } else {
      if (viewPtIdx >= (int) viewPtList.size())
         viewPtIdx--;
      writeViewPtIdx();
      setViewPt();
   }
}


// Renames currently selected viewpoint.

void G4OpenInventorXtExaminerViewer::renameViewPt(char *vpName)
{
   int idx = 0, end, pos;
   int size = viewPtList.size();
   std::string line, newName;
   fileIn.open(fileName.c_str());

   newName = vpName;
   while ((int) newName.size() < MAX_VP_NAME)
      newName += " ";

   getline(fileIn, line);
   pos = fileIn.tellg();
   while (getline(fileIn, line)) {
      end = line.find_last_not_of(' ');
      line = line.substr(0, end + 1);
      if (!strcmp(line.c_str(), curViewPtName)) {
         fileOut.seekp(pos);
         fileOut << newName;
         fileOut.seekp(0, std::ios::end); // Set the file pointer to the end of the file
         break;
      }
      while (line.size())
         getline(fileIn, line);
      pos = fileIn.tellg();
   }

   fileIn.close();
   fileIn.clear();

   while (idx < size) {
      if (!strcmp(viewPtList[idx].viewPtName, curViewPtName)) {
         strcpy(viewPtList[idx].viewPtName, vpName);
         break;
      }
      idx++;
   }
}


// Rewrites entire viewpoint file with sorted viewpoints.

void G4OpenInventorXtExaminerViewer::sortViewPts(std::vector<std::string> sortedViewPts) 
{
   SbVec3f axis;
   float x, y, z, angle;
   int sortIdx = 0, unsortIdx = 0;

   if (fileOut.is_open())
      fileOut.close();

   fileOut.open(fileName.c_str()); // Erase current viewpoint file

   writeViewPtIdx();

   int size = sortedViewPts.size();
   while (sortIdx < size) {
      while (strcmp(sortedViewPts[sortIdx].c_str(),
                    viewPtList[unsortIdx].viewPtName))
         unsortIdx++;

      std::string vpName = viewPtList[unsortIdx].viewPtName;

      while ((int) vpName.size() < MAX_VP_NAME)
         vpName += " ";
      fileOut << vpName << std::endl;
      viewPtList[unsortIdx].position.getValue(x, y, z);
      fileOut << x << " " << y << " " << z << std::endl;

      // Reusing x, y and z for storing the axis
      viewPtList[unsortIdx].orientation.getValue(axis, angle);
      axis.getValue(x, y, z);
      fileOut << x << " " << y << " " << z << " " << angle << std::endl;

      fileOut << viewPtList[unsortIdx].camType << " "
              << viewPtList[unsortIdx].height << std::endl;
      fileOut << viewPtList[unsortIdx].focalDistance << " ";

      fileOut << viewPtList[unsortIdx].nearDistance << " ";

      fileOut << viewPtList[unsortIdx].farDistance << std::endl;

      fileOut << viewPtList[unsortIdx].viewportMapping << " ";
      fileOut << viewPtList[unsortIdx].aspectRatio << "\n" << std::endl;
      fileOut.flush();

      unsortIdx = 0;
      sortIdx++;
   }
}


//  Loads view point data from a file into a vector.

bool G4OpenInventorXtExaminerViewer::loadViewPts() 
{
   bool error = false;
   viewPtData tmp;
   std::string token;
   SbVec3f axis;
   SbRotation orient;
   float x, y, z, angle;

   // Gets the last view point accessed, stored in the first line of the data file.
   fileIn >> token;
   parseString<int>(viewPtIdx, token, error);
   getline(fileIn, token); // Remove "\n"
   // Converts data from string type into necessary types
   while (getline(fileIn, token)) {

      int end = token.find_last_not_of(' '); // Remove padded spaces
      token = token.substr(0, end + 1);

      char *vpName = new char[token.size() + 1];
      strcpy(vpName, token.c_str());
      tmp.viewPtName = vpName;
      fileIn >> token;

      parseString<float>(x, token, error);
      fileIn >> token;
      parseString<float>(y, token, error);
      fileIn >> token;
      parseString<float>(z, token, error);
      fileIn >> token;
      tmp.position = axis.setValue(x, y, z);

      parseString<float>(x, token, error);
      fileIn >> token;
      parseString<float>(y, token, error);
      fileIn >> token;
      parseString<float>(z, token, error);
      fileIn >> token;
      parseString<float>(angle, token, error);
      fileIn >> token;
      orient.setValue(axis.setValue(x, y, z), angle);
      tmp.orientation = orient.getValue();

      int camType;
      parseString<int>(camType, token, error);
      fileIn >> token;
      tmp.camType = (CameraType) camType;

      parseString<float>(tmp.height, token, error);
      fileIn >> token;
      parseString<float>(tmp.focalDistance, token, error);
      fileIn >> token;
      parseString<float>(tmp.nearDistance, token, error);
      fileIn >> token;
      parseString<float>(tmp.farDistance, token, error);
      fileIn >> token;
      parseString<int>(tmp.viewportMapping, token, error);
      fileIn >> token;
      parseString<float>(tmp.aspectRatio, token, error);

      getline(fileIn, token); // To remove "\n" characters
      getline(fileIn, token);

      if (error) {
         viewPtIdx = 0;
         viewPtList.clear();
         return false;
      }
      viewPtList.push_back(tmp);
   }

   return true;
}


// Converts a string type word into a float type.

template<class T> 
void G4OpenInventorXtExaminerViewer::parseString(T &t, const std::string &s,
                                                 bool &error) 
{
   std::istringstream str(s);
   if ((str >> t).fail())
      error = true;
}


// Generic fileSelectionDialog creation.

void G4OpenInventorXtExaminerViewer::popUpFileSelDialog(Widget &dialog,
                                                std::string dialogName,
                                                std::string buttonLabel,
                                                XtCallbackProc cbOK)
{
   int n;
   Arg args[3];
   Widget parent, scrollWidget;
   parent = SoXt::getShellWidget(getParentWidget());

   if (dialog == NULL) {

      // Change the 'OK' button to whatever buttonLabel contains
      XmString str = XmStringCreateLocalized((char *) buttonLabel.c_str());

      n = 0;
      XtSetArg(args[n], XmNokLabelString, str);		n++;
      XtSetArg(args[n], XmNresizePolicy, XmRESIZE_NONE);		n++;

      dialog = XmCreateFileSelectionDialog(parent,
                                           (char *) dialogName.c_str(), args, n);

      XtAddCallback(dialog, XmNokCallback, cbOK, this);
      XtAddCallback(dialog, XmNcancelCallback, cancelFileSelDialogCB, this);

      // Adding scrolling functionality to the widget
      scrollWidget = XmFileSelectionBoxGetChild(dialog, XmDIALOG_DIR_LIST);
      if (scrollWidget)
         xmAddMouseEventHandler(scrollWidget);
      scrollWidget = XmFileSelectionBoxGetChild(dialog, XmDIALOG_LIST);
      if (scrollWidget)
         xmAddMouseEventHandler(scrollWidget);

      XtUnmanageChild(XmSelectionBoxGetChild(dialog, XmDIALOG_HELP_BUTTON));
      XmStringFree(str);
   }
   XtManageChild(dialog);
}


// Generic fileSelectionDialog cancelation.

void G4OpenInventorXtExaminerViewer::cancelFileSelDialogCB(Widget w,
                                                           XtPointer,
                                                           XtPointer)
{
   XtUnmanageChild(w);
}


// Displays a file selection dialog that allows to open a new viewpoint file.

void G4OpenInventorXtExaminerViewer::openViewPtFileCB(Widget,
                                       XtPointer client_data,
                                       XtPointer)
{
   G4OpenInventorXtExaminerViewer * This =
      (G4OpenInventorXtExaminerViewer *) client_data;
   This->popUpFileSelDialog(This->openFileDialog, "Open File", "Load",
                            viewPtFileSelectedCB);
}


void G4OpenInventorXtExaminerViewer::viewPtFileSelectedCB(Widget w, 
                                             XtPointer client_data,
                                             XtPointer call_data)
{
   char *file = NULL;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   XmFileSelectionBoxCallbackStruct *cbs =
      (XmFileSelectionBoxCallbackStruct *) call_data;

   // Get the file
   if (cbs) {
      if (!(file = (char *) XmStringUnparse(cbs->value,
                                            XmFONTLIST_DEFAULT_TAG, XmCHARSET_TEXT, XmCHARSET_TEXT, NULL, 0,
                                            XmOUTPUT_ALL))) {
         SoDebugError::post("G4OpenInventorXtExaminerViewer::fileSelectedCB",
                            "Internal error during file opening");
         return;
      }

      This->fileIn.open(file);
      if (!This->fileIn.fail()) {
         // Opens a file without erasing it
         This->cleanUpAfterPrevFile();
         if (!This->loadViewPts()) {
            String dialogName = (char *) "Error Loading File";
            std::string msg = "Wrong or corrupted input file.";
            This->warningMsgDialog(msg, dialogName, NULL);
         } else {
            This->fileName = file;
            This->fileOut.open(This->fileName.c_str(), std::ios::in);
            This->fileOut.seekp(0, std::ios::end);

            if (!This->listsDialog)
               constructListsDialog(w, This, NULL); // Pop up listsDialog
            else
               This->addViewPoints();

            std::string newDialogName = This->fileName.substr(
                                                              This->fileName.rfind('/') + 1);
            XtVaSetValues(This->myShellDialog, XmNtitle,
                          (char *) newDialogName.c_str(), NULL);

            if (This->viewPtList.size()) {
               This->setViewPt();
               XmTextSetString(This->viewPtSelection, NULL);
               XtSetSensitive(This->nextViewPtButton, True);
               XtSetSensitive(This->prevViewPtButton, True);
            } else {
               XtSetSensitive(This->nextViewPtButton, False);
               XtSetSensitive(This->prevViewPtButton, False);
            }

            XtUnmanageChild(w);
         }

         This->fileIn.close();
      } else {
         String dialogName = (char *) "Nonexistent File";
         std::string msg = "Unable to open file.";
         This->warningMsgDialog(msg, dialogName, NULL);
      }
   }

   This->fileIn.clear();
   XtFree(file);
}


// Adds bookmarks to listsDialog.

void G4OpenInventorXtExaminerViewer::addViewPoints()
{
   int size = viewPtList.size();
   if (!size)
      return;

   XmString *viewPts;

   viewPts = (XmString *) XtMalloc(size * sizeof(XmString));
   for (int i = 0; i < size; i++)
      viewPts[i] = XmStringCreateLocalized(viewPtList[i].viewPtName);

   XmListAddItemsUnselected(myViewPtList, viewPts, size, 1);

   if (viewPts != NULL) {
      for (int i = 0; i < size; i++)
         XmStringFree(viewPts[i]);
      XtFree((char *) viewPts);
   }
}


// Called before loading a new viewpoint file. 
// Resets member fields to default values.

void G4OpenInventorXtExaminerViewer::cleanUpAfterPrevFile()
{
   viewPtIdx = -1;
   viewPtList.clear();
   setSuperimpositionEnabled(superimposition, FALSE);
   scheduleRedraw();
   currentState = GENERAL;
   if (fileOut.is_open())
      fileOut.close();
   if (listsDialog) // Clear viewpoints
      XmListDeleteAllItems(myViewPtList);
}


// Generic function for displaying a warning dialog.

void G4OpenInventorXtExaminerViewer::warningMsgDialog(std::string msg,
                                                      String dialogName,
                                                      XtCallbackProc cb)
{
   Arg args[5];
   unsigned int n;
   XmString warningMsg;

   warningMsg = XmStringCreateLocalized((char *)msg.c_str());

   n = 0;
   XtSetArg(args[n], XmNmessageString, warningMsg);	n++;
   Widget warningDialog = XmCreateWarningDialog(getParentWidget(), dialogName, args, n);
   if (cb)
      XtAddCallback(warningDialog, XmNokCallback, cb, this);

   XmStringFree(warningMsg);

   XtVaSetValues (warningDialog, XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, NULL);
   XtUnmanageChild(XtNameToWidget(warningDialog, "Help"));
   XtUnmanageChild(XtNameToWidget(warningDialog, "Cancel"));

   XtManageChild(warningDialog);
}


void G4OpenInventorXtExaminerViewer::newViewPtFileCB(Widget,
                                                     XtPointer client_data,
                                                     XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = 
      (G4OpenInventorXtExaminerViewer *) client_data;
   This->popUpFileSelDialog(This->newFileDialog, "New File", "Save",
                            createNewVPFileCB);
}


void G4OpenInventorXtExaminerViewer::createNewVPFileCB(Widget w, 
                                                       XtPointer client_data,
                                                       XtPointer call_data)
{
   char *file;
   std::string fName;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   XmFileSelectionBoxCallbackStruct *cbs =
      (XmFileSelectionBoxCallbackStruct *) call_data;

   // Get the file
   if (cbs) {
      if (!(file = (char *) XmStringUnparse(cbs->value,
                                            XmFONTLIST_DEFAULT_TAG, XmCHARSET_TEXT, XmCHARSET_TEXT, NULL, 0,
                                            XmOUTPUT_ALL))) {
         SoDebugError::post("G4OpenInventorXtExaminerViewer::createNewVPFileCB",
                            "Internal error during file opening");
         return;
      }

      This->fileName = file;
      fName = This->fileName.substr(This->fileName.rfind('/') + 1); // Extracts just the name of the file
      This->fileIn.open(file);
      if (This->fileIn.fail()) { // Filename does not exist
         This->cleanUpAfterPrevFile();
         This->fileOut.open(file); // Creates a new empty file
         XtSetSensitive(This->nextViewPtButton, False);
         XtSetSensitive(This->prevViewPtButton, False);
         if (This->listsDialog)
            closeListsDialogCB(w, This, NULL);
         constructListsDialog(w, This, NULL);
         XtUnmanageChild(w);
         if (This->returnToSaveVP) {
            This->returnToSaveVP = false;
            saveViewPtCB(NULL, This, NULL);
         }
      } else { // Filename already exists
         String dialogName = (char *) "Existing File";
         std::string msg = "'" + fName + "' already exists. Do you want to overwrite it?";
         This->warningMsgDialog(msg, dialogName, overwriteFileCB);
         This->fileIn.close();
      }
      This->fileIn.clear();
      XtFree(file);
   }
}


void G4OpenInventorXtExaminerViewer::overwriteFileCB(Widget, 
                                                     XtPointer client_data,
                                                     XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   This->cleanUpAfterPrevFile();
   XtSetSensitive(This->nextViewPtButton, False);
   XtSetSensitive(This->prevViewPtButton, False);

   XtUnmanageChild(This->newFileDialog);

   This->fileOut.open(This->fileName.c_str());

   if (This->returnToSaveVP) {
      This->returnToSaveVP = false;
      saveViewPtCB(NULL, This, NULL);
   }
}


void G4OpenInventorXtExaminerViewer::loadRefCoordsDialogCB(Widget,
                                            XtPointer client_data,
                                            XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = 
      (G4OpenInventorXtExaminerViewer *)client_data;
   This->popUpFileSelDialog(This->loadRefCoordsDialog, "Load Ref Coords",
                            "Load", loadRefCoordsCB);
}


void G4OpenInventorXtExaminerViewer::loadRefCoordsCB(Widget w, 
                                                     XtPointer client_data,
                                                     XtPointer call_data)
{
   char *file = NULL;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *)client_data;
   XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *)call_data;

   // Get the file
   if(cbs) {

      file = (char *)XmStringUnparse(cbs->value, XmFONTLIST_DEFAULT_TAG,
                                     XmCHARSET_TEXT, XmCHARSET_TEXT,
                                     NULL, 0, XmOUTPUT_ALL);

      std::ifstream ifs(file);
      if(ifs.is_open()){
         This->refParticleTrajectory.clear();
         float x,y,z;
         while(ifs >> x >> y >> z){
            This->refParticleTrajectory.push_back(SbVec3f(x,y,z));
         }
         ifs.close();
         XtUnmanageChild(w);
      }
      else{
         String dialogName = (char *) "Problem reading file";
         std::string msg = "Problem reading file";
         This->warningMsgDialog(msg, dialogName, NULL);
         return;

      }
   }

   return;
}


void G4OpenInventorXtExaminerViewer::saveRefCoordsDialogCB(Widget,
                                            XtPointer client_data,
                                            XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (!This->refParticleTrajectory.size()) {
      String dialogName = (char *) "No Reference Trajectory";
      std::string msg = "You need to start a run or load a reference trajectory from a file";
      This->warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   int n;
   Arg args[3];
   Widget parent, scrollWidget;
   parent = SoXt::getShellWidget(This->getParentWidget());

   if (This->saveRefCoordsDialog == NULL) {

      // Change the 'OK' button to whatever buttonLabel contains
      XmString str = XmStringCreateLocalized((char *)"Save");

      n = 0;
      XtSetArg(args[n], XmNokLabelString, str);		n++;
      XtSetArg(args[n], XmNresizePolicy, XmRESIZE_NONE);		n++;

      This->saveRefCoordsDialog = XmCreateFileSelectionDialog(parent,(char *)"Save Ref Coords", args, n);

      XtAddCallback(This->saveRefCoordsDialog, XmNokCallback, saveRefCoordsCB, This);
      XtAddCallback(This->saveRefCoordsDialog, XmNcancelCallback, cancelFileSelDialogCB, This);

      // Adding scrolling functionality to the widget
      scrollWidget = XmFileSelectionBoxGetChild(This->saveRefCoordsDialog, XmDIALOG_DIR_LIST);
      if (scrollWidget)
         xmAddMouseEventHandler(scrollWidget);
      scrollWidget = XmFileSelectionBoxGetChild(This->saveRefCoordsDialog, XmDIALOG_LIST);
      if (scrollWidget)
         xmAddMouseEventHandler(scrollWidget);

      XtUnmanageChild(XmSelectionBoxGetChild(This->saveRefCoordsDialog, XmDIALOG_HELP_BUTTON));
      XmStringFree(str);
   }

   //TODO: Auto name?

   XtManageChild(This->saveRefCoordsDialog);

}


void G4OpenInventorXtExaminerViewer::saveRefCoordsCB(Widget w, 
                                                     XtPointer client_data,
                                                     XtPointer call_data)
{
   char *file;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   XmFileSelectionBoxCallbackStruct *cbs =
      (XmFileSelectionBoxCallbackStruct *) call_data;

   // Get the file
   if (cbs) {

      file = (char *)XmStringUnparse(cbs->value, XmFONTLIST_DEFAULT_TAG,
                                     XmCHARSET_TEXT, XmCHARSET_TEXT,
                                     NULL, 0, XmOUTPUT_ALL);

      std::ifstream ifile(file);
      if (ifile) {
         //File already exists

         Arg args[4];
         Widget parent = This->getParentWidget();	//gets the dialogshell of the ExaminerViewer widget
         Widget confirmOverwriteDialog;
         XmString msg;

         confirmOverwriteDialog = XmCreateQuestionDialog (parent, (char *)"Confirm overwrite", args, 0);
         msg = XmStringCreateLocalized ((char *)"File exists. Overwrite?");
         XtVaSetValues (confirmOverwriteDialog, XmNmessageString, msg, NULL);

         // If users presses OK, we want to return to this function and
         // save the file.  For that to work, pass it the current widget
         // to be able to grab the filename.
         XtVaSetValues (confirmOverwriteDialog, XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, NULL);
         XtAddCallback (confirmOverwriteDialog, XmNokCallback, saveRefCoordsOverWriteCB, client_data);
         XtAddCallback (confirmOverwriteDialog, XmNcancelCallback, saveRefCoordsOverWriteCB, client_data);

         XmStringFree (msg);

         //The confirmOverwriteDialog will need this
         This->saveRefCoordsFileName = file;
         This->saveRefCoordsWidget = w;

         XtUnmanageChild(XtNameToWidget(confirmOverwriteDialog, "Help"));
         XtManageChild(confirmOverwriteDialog);

         return;
      }
      else{

         std::ofstream ofs(file);
         if(ofs.is_open()){
            float x,y,z;
            for(unsigned int i=0; i < This->refParticleTrajectory.size(); ++i){
               This->refParticleTrajectory[i].getValue(x,y,z);
               ofs << x << " " << y << " " << z << "\n";
            }
            ofs.close();
            XtUnmanageChild(w);
         }
         else{
            String dialogName = (char *) "Error opening file";
            std::string msg = "There was a problem trying to open the file '";
            msg += This->saveRefCoordsFileName;
            msg += "'";

            This->warningMsgDialog(msg, dialogName, NULL);
         }
      }
   }

   return;
}


void G4OpenInventorXtExaminerViewer::saveRefCoordsOverWriteCB(Widget w, 
                                                 XtPointer client_data,
                                                 XtPointer call_data)
{
   XmAnyCallbackStruct *cbs = (XmAnyCallbackStruct *) call_data;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   switch (cbs->reason) {
   case XmCR_OK:
      {
         // Overwrite confirmed, save file and dismiss both
         // dialogs (file dialog and overwrite confirmation dialog)
         std::ofstream ofs(This->saveRefCoordsFileName.c_str());
         if(ofs.is_open()){
            float x,y,z;
            for(unsigned int i=0; i < This->refParticleTrajectory.size(); ++i){
               This->refParticleTrajectory[i].getValue(x,y,z);
               ofs << x << " " << y << " " << z << "\n";
            }
            ofs.close();
            XtUnmanageChild(w);
            XtUnmanageChild(This->saveRefCoordsWidget);
         }
         else{
            String dialogName = (char *) "Error opening file";
            std::string msg = "There was a problem trying to open the file '";
            msg += This->saveRefCoordsFileName;
            msg += "'";

            This->warningMsgDialog(msg, dialogName, NULL);
         }
         break;
      }
   case XmCR_CANCEL:
      {
         // Overwrite refused, dismiss overwrite confirmation
         // dialog and return to file dialog

         // Give focus to the text field instead of the OK button
         XmProcessTraversal(XtNameToWidget(This->saveRefCoordsWidget, "Text"), XmTRAVERSE_CURRENT);

         XtUnmanageChild(w);
         This->saveRefCoordsFileName.clear();
         This->saveRefCoordsWidget = NULL;
         break;
      }
   default:
      return;
   }
}


void G4OpenInventorXtExaminerViewer::loadSceneGraphDialogCB(Widget,
                                             XtPointer client_data,
                                             XtPointer)
{
   G4OpenInventorXtExaminerViewer * This =
      (G4OpenInventorXtExaminerViewer *)client_data;
   This->popUpFileSelDialog(This->loadSceneGraphDialog, "Load Scene Graph",
                            "Load", loadSceneGraphCB);
   return;
}


void G4OpenInventorXtExaminerViewer::loadSceneGraphCB(Widget w, 
                                                      XtPointer client_data,
                                                      XtPointer call_data)
{
   char *file = NULL;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *)client_data;
   XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *)call_data;

   if(cbs) {

      file = (char *)XmStringUnparse(cbs->value, XmFONTLIST_DEFAULT_TAG,
                                     XmCHARSET_TEXT, XmCHARSET_TEXT,
                                     NULL, 0, XmOUTPUT_ALL);

      SoInput sceneInput;
      if (!sceneInput.openFile(file)) {
         String dialogName = (char *) "Problem opening file";
         std::string msg = "Cannot open file ";
         msg += file;
         This->warningMsgDialog(msg, dialogName, NULL);

         sceneInput.closeFile();
         XtUnmanageChild(w);
      }
      // Read the whole file into the database
      This->newSceneGraph = SoDB::readAll(&sceneInput);
      if (This->newSceneGraph == NULL) {
         String dialogName = (char *) "Problem reading file";
         std::string msg = "Problem reading file";
         This->warningMsgDialog(msg, dialogName, NULL);
         return;
      }

      //This->newSceneGraph->ref();
      This->setSceneGraph(This->newSceneGraph);
   }

   return;
}


void G4OpenInventorXtExaminerViewer::saveSceneGraphDialogCB(Widget, 
                                             XtPointer client_data,
                                             XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   int n;
   Arg args[3];
   Widget parent, scrollWidget;
   parent = SoXt::getShellWidget(This->getParentWidget());

   if (This->saveSceneGraphDialog == NULL) {

      // Change the 'OK' button to whatever buttonLabel contains
      XmString str = XmStringCreateLocalized((char *)"Save");

      n = 0;
      XtSetArg(args[n], XmNokLabelString, str);		n++;
      XtSetArg(args[n], XmNresizePolicy, XmRESIZE_NONE);		n++;

      This->saveSceneGraphDialog = XmCreateFileSelectionDialog(parent,(char *)"Save Scene Graph", args, n);

      XtAddCallback(This->saveSceneGraphDialog, XmNokCallback, saveSceneGraphCB, This);
      XtAddCallback(This->saveSceneGraphDialog, XmNcancelCallback, cancelFileSelDialogCB, This);

      // Adding scrolling functionality to the widget
      scrollWidget = XmFileSelectionBoxGetChild(This->saveSceneGraphDialog, XmDIALOG_DIR_LIST);
      if (scrollWidget)
         xmAddMouseEventHandler(scrollWidget);
      scrollWidget = XmFileSelectionBoxGetChild(This->saveSceneGraphDialog, XmDIALOG_LIST);
      if (scrollWidget)
         xmAddMouseEventHandler(scrollWidget);

      XtUnmanageChild(XmSelectionBoxGetChild(This->saveSceneGraphDialog, XmDIALOG_HELP_BUTTON));
      XmStringFree(str);
   }

   //TODO: Auto name?

   XtManageChild(This->saveSceneGraphDialog);

}



void G4OpenInventorXtExaminerViewer::saveSceneGraphCB(Widget w, 
                                                      XtPointer client_data,
                                                      XtPointer call_data)
{
   char *file;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   XmFileSelectionBoxCallbackStruct *cbs =
      (XmFileSelectionBoxCallbackStruct *) call_data;

   if (cbs) {

      file = (char *)XmStringUnparse(cbs->value, XmFONTLIST_DEFAULT_TAG,
                                     XmCHARSET_TEXT, XmCHARSET_TEXT,
                                     NULL, 0, XmOUTPUT_ALL);

      std::ifstream ifile(file);
      if (ifile) {
         //File already exists

         Arg args[4];
         Widget parent = This->getParentWidget();	//gets the dialogshell of the ExaminerViewer widget
         Widget confirmOverwriteDialog;
         XmString msg;

         confirmOverwriteDialog = XmCreateQuestionDialog (parent, (char *)"Confirm overwrite", args, 0);
         msg = XmStringCreateLocalized ((char *)"File exists. Overwrite?");
         XtVaSetValues (confirmOverwriteDialog, XmNmessageString, msg, NULL);

         // If users presses OK, we want to return to this function and
         // save the file.  For that to work, pass it the current widget
         // to be able to grab the filename.
         XtVaSetValues (confirmOverwriteDialog, XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, NULL);
         XtAddCallback (confirmOverwriteDialog, XmNokCallback, saveSceneGraphOverWriteCB, client_data);
         XtAddCallback (confirmOverwriteDialog, XmNcancelCallback, saveSceneGraphOverWriteCB, client_data);

         XmStringFree (msg);

         //The confirmOverwriteDialog will need this
         This->saveScenegraphFileName = file;
         This->saveScenegraphWidget = w;

         XtUnmanageChild(XtNameToWidget(confirmOverwriteDialog, "Help"));
         XtManageChild(confirmOverwriteDialog);

         return;
      }
      else{

         SoWriteAction writeAction;
         SoSeparator *root = (SoSeparator *) (This->getSceneGraph());

         SoOutput * out = writeAction.getOutput();

         if(out->openFile(file)){
            out->setBinary(FALSE);
            writeAction.apply(root);
            out->closeFile();

            XtUnmanageChild(w);
         }
         else{
            String dialogName = (char *) "Error opening file";
            std::string msg = "There was a problem trying to open the file '";
            msg += This->saveScenegraphFileName;
            msg += "'";

            This->warningMsgDialog(msg, dialogName, NULL);
         }

      }
   }

   return;
}



void G4OpenInventorXtExaminerViewer::saveSceneGraphOverWriteCB(Widget w, 
                                                  XtPointer client_data,
                                                  XtPointer call_data)
{
   XmAnyCallbackStruct *cbs = (XmAnyCallbackStruct *) call_data;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   switch (cbs->reason) {
   case XmCR_OK:
      {
         // Overwrite confirmed, save file and dismiss both
         // dialogs (file dialog and overwrite confirmation dialog)
         SoWriteAction writeAction;
         SoSeparator *root = (SoSeparator *) (This->getSceneGraph());

         SoOutput * out = writeAction.getOutput();
         if(out->openFile(This->saveScenegraphFileName.c_str())){
            out->setBinary(FALSE);
            writeAction.apply(root);
            out->closeFile();

            XtUnmanageChild(w);
            XtUnmanageChild(This->saveScenegraphWidget);
            This->saveScenegraphFileName.clear();
            This->saveScenegraphWidget = NULL;
         }
         else{
            String dialogName = (char *) "Error opening file";
            std::string msg = "There was a problem trying to open the file '";
            msg += This->saveScenegraphFileName;
            msg += "'";

            This->warningMsgDialog(msg, dialogName, NULL);
            This->saveScenegraphFileName.clear();
            This->saveScenegraphWidget = NULL;
         }
         break;
      }
   case XmCR_CANCEL:
      {
         // Overwrite refused, dismiss overwrite confirmation
         // dialog and return to file dialog

         // Give focus to the text field instead of the OK button
         XmProcessTraversal(XtNameToWidget(This->saveScenegraphWidget, "Text"), XmTRAVERSE_CURRENT);

         XtUnmanageChild(w);
         This->saveScenegraphFileName.clear();
         This->saveScenegraphWidget = NULL;
         break;
      }
   default:
      return;
   }
}


// Receives the name of the bookmark clicked and searches for it in viewPtList.

void G4OpenInventorXtExaminerViewer::loadBookmarkCB(Widget, 
                                                    XtPointer client_data,
                                                    XtPointer call_data)
{
   char *vpName;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;
   XmListCallbackStruct *cbs = (XmListCallbackStruct *) call_data;

   vpName = (char *) XmStringUnparse(cbs->item, XmFONTLIST_DEFAULT_TAG,
                                     XmCHARSET_TEXT, XmCHARSET_TEXT, NULL, 0, XmOUTPUT_ALL);

   for (int i = 0; i < (int) This->viewPtList.size(); i++) {
      if (!strcmp(This->viewPtList[i].viewPtName, vpName)) {
         This->viewPtIdx = i;
         break;
      }
   }
   XmTextSetString(This->viewPtSelection, vpName);

   This->writeViewPtIdx();
   This->setViewPt();
   XtFree(vpName);
}



void G4OpenInventorXtExaminerViewer::deleteBookmarkCB(Widget, 
                                                      XtPointer client_data,
                                                      XtPointer)
{
   char *vpName;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   vpName = XmTextGetString(This->viewPtSelection);

   XmString vpNameStr = XmStringCreateLocalized(vpName);

   if (XmListItemExists(This->myViewPtList, vpNameStr)) {
      XmListDeleteItem(This->myViewPtList, vpNameStr);
      This->deleteViewPt(vpName);
   }

   XmStringFree(vpNameStr);
   XmTextSetString(This->viewPtSelection, NULL);
   XtFree(vpName);
}


void G4OpenInventorXtExaminerViewer::renameBookmarkCB(Widget, 
                                                      XtPointer client_data,
                                                      XtPointer)
{
   std::string vpNameStr;
   char *vpName;
   int *pos_list, pos_cnt;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   vpName = XmTextGetString(This->viewPtSelection);

   if (!strlen(vpName) || !strcmp(This->curViewPtName, vpName)) {
      XtFree(vpName);
      return;
   }

   vpNameStr = vpName;
   XtFree(vpName);
   int beg = vpNameStr.find_first_not_of(' '); // Remove leading/trailing spaces
   int end = vpNameStr.find_last_not_of(' ');
   vpNameStr = vpNameStr.substr(beg, end - beg + 1);
   const int nVPName = vpNameStr.size() + 1;
   char* vpName1 = new char[nVPName];
   strncpy(vpName1, vpNameStr.c_str(), nVPName);

   int size = This->viewPtList.size();
   for (int i = 0; i < size; i++) {
      if (!strcmp(vpName1, This->viewPtList[i].viewPtName)) {

         String dialogName = (char *) "Existing Viewpoint";
         std::string msg = "'";
         msg += vpName1;
         msg += "' already exists. Choose a different name";

         This->warningMsgDialog(msg, dialogName, NULL);
         delete[] vpName1;
         return;
      }
   }

   XmString vpNameXmStr = XmStringCreateLocalized(vpName1);

   if (XmListGetSelectedPos(This->myViewPtList, &pos_list, &pos_cnt)) {
      XmListReplaceItemsPos(This->myViewPtList, &vpNameXmStr, 1, pos_list[0]);
      This->renameViewPt(vpName1);
      XtFree((char *) pos_list);
   }

   if (This->currentState == VIEWPOINT)
      This->scheduleRedraw();

   XmStringFree(vpNameXmStr);
   delete[] vpName1;
}


void G4OpenInventorXtExaminerViewer::sortBookmarksCB(Widget, 
                                                     XtPointer client_data,
                                                     XtPointer)
{
   int size;
   char *vpName;
   XmString *strList, *newStrList;
   std::vector<std::string> charList;
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (This->viewPtList.size() < 2)
      return;

   // Get current entries from the list
   XtVaGetValues(This->myViewPtList, XmNitemCount, &size, XmNitems, &strList,
                 NULL);

   for (int i = 0; i < size; i++) {
      vpName = (char *) XmStringUnparse(strList[i], XmFONTLIST_DEFAULT_TAG,
                                        XmCHARSET_TEXT, XmCHARSET_TEXT, NULL, 0, XmOUTPUT_ALL);
      charList.push_back(vpName);
      XtFree(vpName);
   }

   std::sort(charList.begin(), charList.end());

   newStrList = (XmString *) XtMalloc(size * sizeof(XmString));
   for (int i = 0; i < size; i++) {
      // viewPtIdx has to be changed to account for a different order in viewPtList
      if (!strcmp(charList[i].c_str(), This->curViewPtName))
         This->viewPtIdx = i;
      const int nVPName = charList[i].size() + 1;
      char *vpName2 = new char[nVPName];
      strncpy(vpName2, charList[i].c_str(), nVPName);
      newStrList[i] = XmStringCreateLocalized(vpName2);
      delete [] vpName2;
   }

   XmListDeleteAllItems(This->myViewPtList);
   XmListAddItemsUnselected(This->myViewPtList, newStrList, size, 1);

   This->sortViewPts(charList);

   if (newStrList != NULL) {
      for (int i = 0; i < size; i++)
         XmStringFree(newStrList[i]);
      XtFree((char *) newStrList);
   }
}


void G4OpenInventorXtExaminerViewer::evenOutRefParticlePts()
{
   if(this->refParticleTrajectory.empty())
      return;

   SbVec3f p1, p2, p3, dirNow, dirNxt, dir, p2_tmp, p_start, p_corner, p_nxt;
   float avgDistBtwPts = 0;
   float totalDistBtwPts = 0;
   std::vector<SbVec3f> newRefParticleTrajectory;
   SbVec3f refPoint;
   int size = refParticleTrajectory.size() - 1;
   int numOfPts = 0;
   for (int i = 0; i < size; i++) {
      p1 = refParticleTrajectory[i];
      p2 = refParticleTrajectory[i + 1];
      if (p1 == p2)
         continue;
      numOfPts++;
      totalDistBtwPts += (p2 - p1).length();
   }

   avgDistBtwPts = totalDistBtwPts / numOfPts;
   float minDistAllowed = 0.75 * avgDistBtwPts;
   //	float maxDistAllowed = 1.25 * avgDistBtwPts; // Pts tend to be close not far

   float x, y, z;
   int i = 0, j = 0;
   while (i < size) {
      p1 = refParticleTrajectory[i];
      p2 = refParticleTrajectory[i + 1];

      refPoint = p1;
      p1.getValue(x, y, z);

      newRefParticleTrajectory.push_back(refPoint);

      j = i;
      while ((p2 - p1).length() < minDistAllowed && j < (size - 1)) {
         j++;

         p1 = refParticleTrajectory[j];
         p2 = refParticleTrajectory[j + 1];
      }
      if (j != i)
         i = j + 1;
      else
         i++;
   }

   refParticleTrajectory.clear();
   refParticleTrajectory = newRefParticleTrajectory;
}


// Called when the viewer is closed; closes all open widgets.

void G4OpenInventorXtExaminerViewer::closeMainWindowCB(Widget, 
                                                       XtPointer client_data,
                                                       XtPointer)
{
   G4OpenInventorXtExaminerViewer * This =
      (G4OpenInventorXtExaminerViewer *) client_data;

   if (This->openFileDialog)
      XtUnmanageChild(This->openFileDialog);

   if (This->newFileDialog)
      XtUnmanageChild(This->newFileDialog);

   if (This->listsDialog)
      closeListsDialogCB(NULL, This, NULL);
}


void G4OpenInventorXtExaminerViewer::saveCurCamera()
{
   SoCamera *cam = getCamera();
   camB4Animation.viewportMapping = cam->viewportMapping.getValue();
   camB4Animation.position = cam->position.getValue();
   camB4Animation.orientation = cam->orientation.getValue();
   camB4Animation.aspectRatio = cam->aspectRatio.getValue();
   camB4Animation.nearDistance = cam->nearDistance.getValue();
   camB4Animation.farDistance = cam->farDistance.getValue();
   camB4Animation.focalDistance = cam->focalDistance.getValue();

   if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
      camB4Animation.height =
         ((SoPerspectiveCamera *) cam)->heightAngle.getValue();
      camB4Animation.camType = PERSPECTIVE;
   } else if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
      camB4Animation.height =
         ((SoOrthographicCamera *) cam)->height.getValue();
      camB4Animation.camType = ORTHOGRAPHIC;
   }
}


void G4OpenInventorXtExaminerViewer::restoreCamera()
{
   SoCamera *cam = getCamera();

   cam->viewportMapping = camB4Animation.viewportMapping;
   cam->position = camB4Animation.position;
   cam->orientation = camB4Animation.orientation;
   cam->aspectRatio = camB4Animation.aspectRatio;
   cam->nearDistance = camB4Animation.nearDistance;
   cam->farDistance = camB4Animation.farDistance;
   cam->focalDistance = camB4Animation.focalDistance;

   if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
      if (camB4Animation.camType == ORTHOGRAPHIC) {
         toggleCameraType();
         cam = getCamera();
         ((SoOrthographicCamera *) cam)->height.setValue(
                                                         camB4Animation.height);
      } else
         ((SoPerspectiveCamera *) cam)->heightAngle.setValue(
                                                             camB4Animation.height);
   } else if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
      if (camB4Animation.camType == PERSPECTIVE) {
         toggleCameraType();
         cam = getCamera();
         ((SoPerspectiveCamera *) cam)->heightAngle.setValue(
                                                             camB4Animation.height);
      } else
         ((SoOrthographicCamera *) cam)->height.setValue(
                                                         camB4Animation.height);
   }
}


void G4OpenInventorXtExaminerViewer::animateSensorRotationCB(void *data, 
                                                             SoSensor *sensor)
{
   SbTime curTime = SbTime::getTimeOfDay();
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) data;
   SoTimerSensor *s = (SoTimerSensor *) sensor;

   float t = float((curTime - s->getBaseTime()).getValue())
      / This->animateBtwPtsPeriod;

   if ((t > 1.0f) || (t + s->getInterval().getValue() > 1.0f))
      t = 1.0f;
   SbBool end = (t == 1.0f);

   if (end) {
      This->animateSensorRotation->unschedule();
      if(This->rotCnt){
         // rotations left
         This->rotateCamera();
      }
      else {
         // rotation over
         This->currentState = This->prevState;
         return;
      }
   }

}


// Called repeatedly during reference particle animation

void G4OpenInventorXtExaminerViewer::animateSensorCB(void *data, 
                                                     SoSensor *sensor)
{
   SbTime curTime = SbTime::getTimeOfDay();
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) data;
   SoCamera *cam = This->getCamera();
   SoTimerSensor *s = (SoTimerSensor *) sensor;

   float t = float((curTime - s->getBaseTime()).getValue())
      / This->animateBtwPtsPeriod;

   if ((t > 1.0f) || (t + s->getInterval().getValue() > 1.0f))
      t = 1.0f;
   SbBool end = (t == 1.0f);

   cam->orientation = SbRotation::slerp(This->camStartOrient, This->camEndOrient, t);
   cam->position = This->camStartPos + (This->camEndPos - This->camStartPos) * t;

   if (end) {
      This->animateSensor->unschedule();

      if (This->currentState == ANIMATION) {
         if (This->refParticleIdx < (int) (This->refParticleTrajectory.size() - 1))
            This->animateRefParticle();
         else {
            This->animateBtwPtsPeriod = MIN_SPEED;
            This->speedStep = START_STEP;
         }
      }
      if (This->currentState == REVERSED_ANIMATION) {
         if (This->refParticleIdx >= 1)
            This->animateRefParticle();
         else {
            This->animateBtwPtsPeriod = MIN_SPEED;
            This->speedStep = START_STEP;
         }
      }
   }
}


void G4OpenInventorXtExaminerViewer::setStartingPtForAnimation()
{
   if (SoXtExaminerViewer::isAnimating())
      stopAnimating();

   SbRotation rot;
   SbVec3f p1, p2, p2_tmp, camUpV, camD, camD_tmp, leftRightAxis;
   float x1, y1, z1, x2, y2, z2;

   if (currentState == ANIMATION) {
      p1 = refParticleTrajectory[refParticleIdx];
      p2 = refParticleTrajectory[++(refParticleIdx)];
   } else if (currentState == REVERSED_ANIMATION) {
      p2 = refParticleTrajectory[refParticleIdx];
      p1 = refParticleTrajectory[--(refParticleIdx)];
   } else if (currentState == PAUSED_ANIMATION) {
      if (refParticleIdx < (int) refParticleTrajectory.size()) {
         p1 = refParticleTrajectory[refParticleIdx];
         p2 = refParticleTrajectory[refParticleIdx + 1];
      } else {
         p1 = refParticleTrajectory[refParticleIdx - 1];
         p2 = refParticleTrajectory[refParticleIdx];
      }
   }
   p1.getValue(x1, y1, z1);
   p2.getValue(x2, y2, z2);

   camD = p2 - p1;
   camD.normalize();

   p2_tmp.setValue(x2, y1, z2);
   camD_tmp = p2_tmp - p1;
   camD_tmp.normalize();

   camUpV.setValue(0, 1, 0);
   rot.setValue(camD_tmp, camD);
   rot.multVec(camUpV, camUpV);

   leftRightAxis = camD.cross(camUpV);

   myCam->position = p1;
   myCam->pointAt(p2, camUpV);

   // Update camera position
   p1 = p1 + (up_down * camUpV) + (left_right * leftRightAxis);
   myCam->position = p1;
   // FWJ Try look-ahead here
   int idx = refParticleIdx + pathLookahead;
   idx = std::min(idx, (int)refParticleTrajectory.size() - 1);
   myCam->pointAt(refParticleTrajectory[idx], camUpV);
   //   myCam->pointAt(refParticleTrajectory[idx], camUpVec);
   myCam->focalDistance = 0.1f;
}


void G4OpenInventorXtExaminerViewer::gotoRefPathStart()
{
   G4OpenInventorXtExaminerViewer::gotoRefPathStartCB(NULL, (void *)this, 
                                                      NULL);
}


void G4OpenInventorXtExaminerViewer::gotoRefPathStartCB(Widget, 
                                                        XtPointer client_data, 
                                                        XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (!This->refParticleTrajectory.size()) {
      String dialogName = (char *) "No Reference Trajectory";
      std::string msg = "You need to start a run or load a reference trajectory from a file";
      This->warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   if (This->currentState == ROTATING)
      return;
   if (This->currentState == ANIMATION || This->currentState == REVERSED_ANIMATION
       || This->currentState == PAUSED_ANIMATION) {
      if (This->animateSensor->isScheduled())
         This->animateSensor->unschedule();
      This->setSuperimpositionEnabled(This->superimposition, FALSE);
      This->maxSpeed = 0.0f;
      This->scheduleRedraw();
   } else {
      This->saveCurCamera();
      This->prevState = This->currentState;
      This->prevRefIdx = This->refParticleIdx;
   }

   if (This->SoXtExaminerViewer::isAnimating())
      This->stopAnimating();

   This->up_down = 0;
   This->left_right = 0;
   This->step = 1;

   This->refParticleIdx = 0;
   This->currentState = BEAMLINE;
   This->setSuperimpositionEnabled(This->superimposition, TRUE);
   This->axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
   This->animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
   This->animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);
   This->scheduleRedraw();

   // FWJ Disabled: this is set in moveCamera()
   // Zoom in
   //   SoCamera *cam = This->getCamera();
   //   cam->focalDistance = 0.1f;

   This->prevParticleDir = SbVec3f(0,0,0);

   //Default zoom
   SbVec3f p1 = This->refParticleTrajectory[0];
   SbVec3f pN = This->refParticleTrajectory[This->refParticleTrajectory.size() - 1];
   This->distance = (pN - p1).length() / 10;

   This->moveCamera(This->distance, true);
}


void G4OpenInventorXtExaminerViewer::invertRefPathCB(Widget,
                                                     XtPointer client_data,
                                                     XtPointer)
{
   G4OpenInventorXtExaminerViewer * This =
      (G4OpenInventorXtExaminerViewer *) client_data;
   This->invertRefPath();
}


void G4OpenInventorXtExaminerViewer::invertRefPath()
{
   std::reverse(this->refParticleTrajectory.begin(),
                this->refParticleTrajectory.end());
   this->setReferencePathZPos();
   this->sortElements();
}


void G4OpenInventorXtExaminerViewer::animateRefParticleCB(Widget, 
                                           XtPointer client_data,
                                           XtPointer)
{
   G4OpenInventorXtExaminerViewer * This = (G4OpenInventorXtExaminerViewer *) client_data;

   if (!This->refParticleTrajectory.size()) {
      This->returnToAnim = true;
      String dialogName = (char *) "No Reference Trajectory";
      std::string msg = "You need to start a run or load a reference trajectory from a file";
      This->warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   if (!This->refParticleTrajectory.size())
      return;

   ///////////////////////////////////////////////////////////////
   This->setSuperimpositionEnabled(This->superimposition, TRUE);
   This->maxSpeed = SPEED_INDICATOR_STEP;
   This->axisSwitch->whichChild.setValue(SO_SWITCH_ALL);
   This->animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_ALL);
   This->animSpeedSwitch->whichChild.setValue(SO_SWITCH_ALL);
   This->scheduleRedraw();
   ///////////////////////////////////////////////////////////////

   SoCamera *cam = This->getCamera();
   //   SbVec3f camDirOld, camDirNew, camDirNew_tmp, camUpVec, P0, P1, P1_tmp;

   if (This->currentState == ANIMATION || This->currentState == REVERSED_ANIMATION
       || This->currentState == ROTATING)
      return;

   if (This->currentState != PAUSED_ANIMATION) {

      This->saveCurCamera();
      This->prevState = This->currentState;
      This->prevRefIdx = This->refParticleIdx;

      if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
         This->toggleCameraType();
         cam = This->getCamera();
      }

      This->refParticleIdx = 0; // Set the camera to the starting point of the animation
      This->animateBtwPtsPeriod = MIN_SPEED;
      This->speedStep = START_STEP;
      This->left_right = This->up_down = 0;

      cam->focalDistance = 0.1f;
      ((SoPerspectiveCamera *) cam)->heightAngle = 0.50f;
   }

   This->currentState = ANIMATION;
   This->setStartingPtForAnimation();

   cam->position = (This->myCam)->position.getValue();
   cam->orientation = (This->myCam)->orientation.getValue();
   This->animateRefParticle(); // Animate the camera
}


void G4OpenInventorXtExaminerViewer::animateRefParticle()
{
   SoCamera *cam = getCamera();

   camStartPos = cam->position.getValue();
   camStartOrient = cam->orientation.getValue();

   if (currentState != BEAMLINE)
      setStartingPtForAnimation();

   camEndPos = myCam->position.getValue();
   camEndOrient = myCam->orientation.getValue();

   if (animateSensor->isScheduled())
      animateSensor->unschedule();

   animateSensor->setBaseTime(SbTime::getTimeOfDay());
   animateSensor->setInterval(SbTime(0.02));

   animateSensor->schedule();
}


void G4OpenInventorXtExaminerViewer::addEscapeCallback(
                         void (*callback)(void *), void * object)
{
   this->escapeCallback = callback;
   this->examinerObject = object;
}


void G4OpenInventorXtExaminerViewer::sceneChangeCB(void *userData, SoSensor *)
{
   G4OpenInventorXtExaminerViewer* This =
      (G4OpenInventorXtExaminerViewer*)userData;
   if(This->newEvents){
      This->findAndSetRefPath();
      This->newEvents = false;
   }
}


HookEventProcState::HookEventProcState(G4OpenInventorXtExaminerViewer* vwr)
{
   this->viewer = vwr;
}


HookEventProcState::~HookEventProcState()
{;}


G4bool HookEventProcState::Notify(G4ApplicationState requiredState)
{
   if(requiredState == G4State_EventProc){
      this->viewer->newEvents = true;
   }
   return true;
}
