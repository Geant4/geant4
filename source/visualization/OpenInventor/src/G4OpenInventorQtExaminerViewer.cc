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

#include "G4OpenInventorQtExaminerViewer.hh"

#include "ui_OIQtListsDialog.h"

#include "saveViewPt.h"
#include "pickext.h"
#include "pickref.h"
#include "wireframe.h"

#include <algorithm> // For using sort on a vector

#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4UIQt.hh"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/SoQtCursor.h>
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
#include <Inventor/projectors/SbPlaneProjector.h>

#include <Inventor/sensors/SoTimerSensor.h>   // Animation
#include <Inventor/sensors/SoNodeSensor.h>    // Detect start of run

#include "Geant4_SoPolyhedron.h"
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

#include <QMenuBar>
#include <QPushButton>
#include <QRadioButton>
#include <QToolButton>
#include <QListWidget>
#include <QListWidgetItem>
#include <QInputDialog>
#include <QMessageBox>
#include <QFileDialog>
#include <QStyle>
#include <QCommonStyle>
//#include <QMainWindow>

#ifndef G4GMAKE
#include "moc_G4OpenInventorQtExaminerViewer.cpp"
#endif

#define G4warn G4cout

G4OpenInventorQtExaminerViewer* G4OpenInventorQtExaminerViewer::viewer = 0;

#define MIN_SPEED  2.1        // Lower number means faster
#define START_STEP 0.3
#define SPEED_INDICATOR_STEP 0.045
#define MAX_SPEED_INDICATOR  0.81
// Number of steps 90 degree rotation around an element is split into
#define ROT_CNT 6


// Constructor
G4OpenInventorQtExaminerViewer::
G4OpenInventorQtExaminerViewer(QWidget* parent, const char* name, SbBool embed,
                               SoQtFullViewer::BuildFlag flag,
                               SoQtViewer::Type type)
   : SoQtExaminerViewer(parent, name, embed, flag, type),
     externalQtApp(0), processSoEventCount(0)
{
   // FWJ DEBUG
   //  G4cout << "G4OpenInventorQtExaminerViewer CONSTRUCTOR CALLED" << G4endl;
   //  G4cout << "G4OpenInventorQtExaminerViewer parent=" << parent << G4endl;

   // FWJ THIS DOESN'T WORK APPARENTLY NO MAINWINDOW
   //   QMenuBar* menubar = ((QMainWindow*)parent)->menuBar();

   fName = new QString(name);
   viewer = this;
   construct(TRUE);
}

// Destructor
G4OpenInventorQtExaminerViewer::~G4OpenInventorQtExaminerViewer()
{
   //   if (superimposition != NULL) {
   //      removeSuperimposition(superimposition);
   //      superimposition->unref();
   //      superimposition = NULL;
   //   }
   //   if (animateSensor->isScheduled())
   //      animateSensor->unschedule();
   //   delete animateSensor;
   //   delete sceneChangeSensor;
   //   delete[] curViewPtName;
   //   delete searcher;

   viewer = 0;
}


void G4OpenInventorQtExaminerViewer::construct(const SbBool)
{
   setFeedbackSize(40);

   hookBeamOn = new HookEventProcState(this);
   newEvents = false;

   buildWidget(getParentWidget());

   fileName = "bookmarkFile"; // Default viewpoint file name
   viewPtIdx = -1; // index of the most recent viewpoint in viewPtList vector

   animateSensor = new SoTimerSensor(animateSensorCB, this);
   animateSensorRotation = new SoTimerSensor(animateSensorRotationCB, this);
   animateBtwPtsPeriod = MIN_SPEED;

   currentState = GENERAL;
   myCam = new SoPerspectiveCamera;
   MAX_VP_IDX = 3;
   MAX_VP_NAME = 35; // Max length of a viewpoint name, padded with spaces
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
   //   warningFlag = false; // We come from the warning dialog
   //   myElementList = NULL;
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
      " MaterialBinding ",
      " {",
      "         value OVERALL",
      " }",
      "         OrthographicCamera ",
      " {",
      "         height 1",
      "         nearDistance 0",
      "         farDistance 1",
      " }",
      "         DEF soxt->callback Callback { }",
      "         Separator ",
      " {",
      "         DEF soxt->translation Translation ",
      "         {",
      "                 translation 0 0 0",
      "     }",
      "     DEF soxt->scale Scale ",
      "         {",
      "                 scaleFactor 1 1 1",
      "     }",
      "         DEF soxt->geometry Coordinate3 ",
      "         {",
      "             point ",
      "                 [",
      "                         -0.81   -0.04   0,      -0.81   0             0,",
      "                 -0.81   0.04    0,      0       -0.04   0,",
      "                 0       0       0,  0       0.04        0,",
      "                 0.81    -0.04   0,  0.81        0           0,",
      "                 0.81    0.04    0,",
      "                 0       0.02    0,", // idx 9
      "                 0.81    0.02    0,  0.81        -0.02   0,",
      "                 0       -0.02   0,",
      "                 0       0.01    0,", // idx 13
      "                 0.4     0.01    0,  0.4         -0.01   0,",
      "                 0       -0.01   0",
      "                 ]",
      "         }",
      // current speed indicator (outline)
      "         DEF soxt->animSpeedOutlineSwitch Switch ",
      "         {",
      "                 whichChild -3",
      "                 Material ",
      "                 {",
      "                    emissiveColor 0 0 0",
      "             }",
      "                 IndexedFaceSet ",
      "                 {",
      "                 coordIndex ",
      "                         [",
      "                                 12, 11, 10, 9, -1",
      "                         ]",
      "         }",
      "                  }",
      // the coordinate system
      "         DEF soxt->axisSwitch Switch ",
      "         {",
      "                 whichChild -3",
      "                 BaseColor ",
      "                 {",
      "                     rgb 1 1 1",
      "                 }",
      "                 IndexedLineSet ",
      "                 {",
      "                         coordIndex ",
      "                         [",
      "                                 0, 2, -1,",
      "                                 3, 5, -1,",
      "                                 6, 8, -1,",
      "                                 1, 7, -1",
      "                         ]",
      "                     }",
      "                 }",
      // current speed indicator
      "         DEF soxt->animSpeedSwitch Switch ",
      "         {",
      "                     whichChild -3",
      "                 Material ",
      "                 {",
      "                 emissiveColor 0 1 0",
      "         }",
      "                 IndexedFaceSet ",
      "                 {",
      "                 coordIndex ",
      "                         [",
      "                                 16, 15, 14, 13, -1",
      "                         ]",
      "                 }",
      "         }",
      "         }",
      // For displaying either z position (during animation) or current viewpoint name
      " DEF soxt->curInfoSwitch Switch ",
      " {",
      "         whichChild -3",
      "         DEF soxt->curInfoTrans Translation ",
      "         {",
      "                 translation 0 0 0    ",
      //      "                 translation 10 20 30    ",
      "         }",
      "         DEF soxt->curInfoFont Font ",
      "         {",
      "                 name defaultFont:Bold",
      "                 size 16",
      "                 }",
      "         DEF soxt->curInfoText Text2 ",
      "         {",
      "                 string Hello",
      "     }",
      " }",
      // Need to use different fields for mouseover
      // because newlines are ignored when the scene is rendered
      " Separator ",
      " {",
      "         DEF soxt->mouseOverTransLogName Translation ",
      "         {",
      "                 translation 0 0 0    ",
      "         }",
      "         DEF soxt->mouseOverFontLogName Font ",
      "         {",
      "                 name defaultFont:Bold",
      "                 size 16",
      "                 }",
      "         DEF soxt->mouseOverTextLogName Text2 { } ",
      " }",
      " Separator ",
      " {",
      "         DEF soxt->mouseOverTransSolid Translation ",
      "         {",
      "                 translation 0 0 0    ",
      "         }",
      "         DEF soxt->mouseOverFontSolid Font ",
      "         {",
      "                 name defaultFont:Bold",
      "                 size 16",
      "                 }",
      "         DEF soxt->mouseOverTextSolid Text2 { } ",
      " }",
      " Separator ",
      " {",
      "         DEF soxt->mouseOverTransMaterial Translation ",
      "         {",
      "                 translation 0 0 0    ",
      "         }",
      "         DEF soxt->mouseOverFontMaterial Font ",
      "         {",
      "                 name defaultFont:Bold",
      "                 size 16",
      "                 }",
      "         DEF soxt->mouseOverTextMaterial Text2 { } ",
      " }",
      " Separator ",
      " {",
      "         DEF soxt->mouseOverTransZPos Translation ",
      "         {",
      "                 translation 0 0 0    ",
      "         }",
      "         DEF soxt->mouseOverFontZPos Font ",
      "         {",
      "                 name defaultFont:Bold",
      "                 size 16",
      "                 }",
      "         DEF soxt->mouseOverTextZPos Text2 { } ",
      " }",
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
   SbBool ok = SoDB::read(input, superimposition);
   (void)ok;   // FWJ added to avoid compiler warning
   assert(ok);
   delete input;
   delete[] buf;
   superimposition->ref();

   sscale = (SoScale *) getSuperimpositionNode(superimposition, "soxt->scale");
   stranslation = (SoTranslation *) getSuperimpositionNode(superimposition, "soxt->translation");
   sgeometry = (SoCoordinate3 *) getSuperimpositionNode(superimposition, "soxt->geometry");
   axisSwitch = (SoSwitch *) getSuperimpositionNode(superimposition, "soxt->axisSwitch");
   animSpeedOutlineSwitch = (SoSwitch *) getSuperimpositionNode(superimposition, "soxt->animSpeedOutlineSwitch");
   animSpeedSwitch = (SoSwitch *) getSuperimpositionNode(superimposition, "soxt->animSpeedSwitch");
   curInfoSwitch = (SoSwitch *) getSuperimpositionNode(superimposition, "soxt->curInfoSwitch");
   curInfoTrans = (SoTranslation *) getSuperimpositionNode(superimposition, "soxt->curInfoTrans");
   curInfoFont = (SoFont *) getSuperimpositionNode(superimposition, "soxt->curInfoFont");
   curInfoText = (SoText2 *) getSuperimpositionNode(superimposition, "soxt->curInfoText");
   mouseOverTransLogName = (SoTranslation*)getSuperimpositionNode(superimposition, "soxt->mouseOverTransLogName");
   mouseOverFontLogName = (SoFont *) getSuperimpositionNode(superimposition, "soxt->mouseOverFontLogName");
   mouseOverTextLogName = (SoText2 *) getSuperimpositionNode(superimposition, "soxt->mouseOverTextLogName");
   mouseOverTransSolid = (SoTranslation *) getSuperimpositionNode(superimposition, "soxt->mouseOverTransSolid");
   mouseOverFontSolid = (SoFont *) getSuperimpositionNode(superimposition, "soxt->mouseOverFontSolid");
   mouseOverTextSolid = (SoText2 *) getSuperimpositionNode(superimposition, "soxt->mouseOverTextSolid");
   mouseOverTransMaterial = (SoTranslation*)getSuperimpositionNode(superimposition, "soxt->mouseOverTransMaterial");
   mouseOverFontMaterial = (SoFont *) getSuperimpositionNode(superimposition, "soxt->mouseOverFontMaterial");
   mouseOverTextMaterial = (SoText2 *) getSuperimpositionNode(superimposition, "soxt->mouseOverTextMaterial");
   mouseOverTransZPos = (SoTranslation *) getSuperimpositionNode(superimposition, "soxt->mouseOverTransZPos");
   mouseOverFontZPos = (SoFont *) getSuperimpositionNode(superimposition, "soxt->mouseOverFontZPos");
   mouseOverTextZPos = (SoText2 *) getSuperimpositionNode(superimposition, "soxt->mouseOverTextZPos");

   SoCallback * cb = (SoCallback *) getSuperimpositionNode(superimposition, "soxt->callback");
   cb->setCallback(superimpositionCB, this);

   addSuperimposition(superimposition);
   setSuperimpositionEnabled(superimposition, FALSE);
   axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
   animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
   animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);

   /////////////////////\SUPERIMPOSED SCENE///////////////////////////////////

}


// Adds a menu bar and menu items to the viewer.
void G4OpenInventorQtExaminerViewer::buildWidget(QWidget* parent)
{
   if (!parent)
      SoDebugError::post("G4OpenInventorQtExaminerViewer::buildWidget",
                         "Error: Parent is null.");

   // Common font for (almost) all widgets
   font = new QFont;
   font->setPointSize(12);
   // This font setting does not propagate to added child widgets - Why?
   parent->setFont(*font);
   // This propagates everywhere but would affect UIQt!
   //   QApplication::setFont(*font);

// MENU BAR

   menubar = new QMenuBar(getRenderAreaWidget());
   // FWJ DEBUG
   //   G4cout << "G4OpenInventorQtExaminerViewer: GOT A menubar=" <<
   //      menubar << G4endl; 
   
   filemenu = new QMenu("File");
   menubar->addMenu(filemenu); 

   FileOpenBookmark = new QAction("Open Bookmark File", this);
   FileOpenBookmark->setFont(*font);
   connect(FileOpenBookmark, SIGNAL(triggered()), this,
           SLOT(FileOpenBookmarkCB()));
   filemenu->addAction(FileOpenBookmark);

   FileNewBookmark = new QAction("New Bookmark File", this);
   FileNewBookmark->setFont(*font);
   connect(FileNewBookmark, SIGNAL(triggered()), this,
           SLOT(FileNewBookmarkCB()));
   filemenu->addAction(FileNewBookmark);

   FileLoadRefPath = new QAction("Load Reference Path", this);
   FileLoadRefPath->setFont(*font);
   connect(FileLoadRefPath, SIGNAL(triggered()), this,
           SLOT(FileLoadRefPathCB()));
   filemenu->addAction(FileLoadRefPath);

   FileSaveRefPath = new QAction("Save Reference Path", this);
   FileSaveRefPath->setFont(*font);
   connect(FileSaveRefPath, SIGNAL(triggered()), this,
           SLOT(FileSaveRefPathCB()));
   filemenu->addAction(FileSaveRefPath);

   FileLoadSceneGraph = new QAction("Load scene graph", this);
   FileLoadSceneGraph->setFont(*font);
   connect(FileLoadSceneGraph, SIGNAL(triggered()), this,
           SLOT(FileLoadSceneGraphCB()));
   filemenu->addAction(FileLoadSceneGraph);

   FileSaveSceneGraph = new QAction("Save scene graph", this);
   FileSaveSceneGraph->setFont(*font);
   connect(FileSaveSceneGraph, SIGNAL(triggered()), this,
           SLOT(FileSaveSceneGraphCB()));
   filemenu->addAction(FileSaveSceneGraph);

   // Rest of File menu is done in G4OpenInventorQtViewer

   toolsmenu = new QMenu("Tools");
   menubar->addMenu(toolsmenu); 

   ToolsAnimateRefParticle = new QAction("Fly on Ref Path", this);
   ToolsAnimateRefParticle->setFont(*font);
   connect(ToolsAnimateRefParticle, SIGNAL(triggered()), this,
           SLOT(ToolsAnimateRefParticleCB()));
   toolsmenu->addAction(ToolsAnimateRefParticle);

   ToolsRefPathStart = new QAction("Go to start of Ref Path", this);
   ToolsRefPathStart->setFont(*font);
   connect(ToolsRefPathStart, SIGNAL(triggered()), this,
           SLOT(ToolsRefPathStartCB()));
   toolsmenu->addAction(ToolsRefPathStart);

   ToolsRefPathInvert = new QAction("Invert Ref Path", this);
   ToolsRefPathInvert->setFont(*font);
   connect(ToolsRefPathInvert, SIGNAL(triggered()), this,
           SLOT(ToolsRefPathInvertCB()));
   toolsmenu->addAction(ToolsRefPathInvert);

   etcmenu = new QMenu("Etc");
   menubar->addMenu(etcmenu); 

   // All Etc menu items are done in G4OpenInventorQtViewer

   helpmenu = new QMenu("Help");
   menubar->addMenu(helpmenu); 

   HelpControls = new QAction("Controls", this);
   HelpControls->setFont(*font);
   connect(HelpControls, SIGNAL(triggered()), this, SLOT(HelpControlsCB()));
   helpmenu->addAction(HelpControls);

   menubar->show();

   //   SoQtExaminerViewer::buildWidget(parent);

   // APP VIEWER BUTTONS have their own box on upper left
   // The built in viewer button list is PRIVATE

   saveViewPtButton = new QPushButton;
   saveViewPtButton->setIcon(QPixmap((const char **)saveViewPt_xpm));
   saveViewPtButton->setIconSize(QSize(24,24));
   saveViewPtButton->setToolTip("Bookmark this view");
   connect(saveViewPtButton, SIGNAL(clicked()), this,
           SLOT(SaveViewPtCB()));
   addAppPushButton(saveViewPtButton);

   nextViewPtButton = new QPushButton;
   nextViewPtButton->setIconSize(QSize(24,24));
   QCommonStyle style;
   nextViewPtButton->setIcon(style.standardIcon(QStyle::SP_ArrowRight));
   nextViewPtButton->setToolTip("Next bookmark");
   connect(nextViewPtButton, SIGNAL(clicked()), this,
           SLOT(NextViewPtCB()));
   addAppPushButton(nextViewPtButton);

   prevViewPtButton = new QPushButton;
   prevViewPtButton->setIconSize(QSize(24,24));
   prevViewPtButton->setIcon(style.standardIcon(QStyle::SP_ArrowLeft));
   prevViewPtButton->setToolTip("Previous bookmark");
   connect(prevViewPtButton, SIGNAL(clicked()), this,
           SLOT(PrevViewPtCB()));
   addAppPushButton(prevViewPtButton);

   abbrOutputButton = new QPushButton;
   abbrOutputButton->setCheckable(true);
   abbrOutputButton->setIconSize(QSize(24,24));
   abbrOutputButton->setIcon(QPixmap((const char **)pickext_xpm));
   abbrOutputButton->setToolTip("Extended picking & readout");
   connect(abbrOutputButton, SIGNAL(toggled(bool)), this,
           SLOT(AbbrOutputCB(bool)));
   addAppPushButton(abbrOutputButton);

   pickRefPathButton = new QPushButton;
   pickRefPathButton->setIconSize(QSize(24,24));
   pickRefPathButton->setIcon(QPixmap((const char **)pickref_xpm));
   pickRefPathButton->setToolTip("Pick ref trajectory");
   connect(pickRefPathButton, SIGNAL(clicked()), this,
           SLOT(PickRefPathCB()));
   addAppPushButton(pickRefPathButton);

   switchWireFrameButton = new QPushButton;
   switchWireFrameButton->setCheckable(true);
   switchWireFrameButton->setIconSize(QSize(24,24));
   switchWireFrameButton->setIcon(QPixmap((const char **)wireframe_xpm));
   switchWireFrameButton->setToolTip("Switch wireframe/solid");
   connect(switchWireFrameButton, SIGNAL(toggled(bool)), this,
           SLOT(SwitchWireFrameCB(bool)));
   addAppPushButton(switchWireFrameButton);

   switchAxesButton = new QPushButton;
   switchAxesButton->setCheckable(true);
   switchAxesButton->setText(QString("A"));
   switchAxesButton->setToolTip("Axes on/off");
   connect(switchAxesButton, SIGNAL(toggled(bool)), this,
           SLOT(SwitchAxesCB(bool)));
   addAppPushButton(switchAxesButton);

   detachButton = new QPushButton;
   detachButton->setIconSize(QSize(24,24));
   detachButton->setIcon(style.standardIcon(QStyle::SP_CommandLink));
   detachButton->setToolTip("Detach viewer window");
   connect(detachButton, SIGNAL(clicked()), this,
           SLOT(DetachCB()));
   // Used for UIQt only so check and add later
   //   addAppPushButton(detachButton);

   // HELP WINDOW

   helpmsgbox = new QMessageBox(getParentWidget());
   helpmsgbox->setWindowTitle("OIQt Controls");
   helpmsgbox->setFont(*font);
   QString messagetxt =
"\nVIEWING mode (Hand cursor):\n\n\
   Left-button + pointer move:  rotate\n\
   Shift+Left-button + pointer move:  pan\n\
   Middle-button + pointer move:  pan\n\
   Ctrl+Shift+Left-button + pointer move:  zoom\n\
   Mouse wheel:  zoom\n\
   Right-button:  popup menu\n\n\
PICKING mode (Arrow cursor):\n\n\
   Click on a volume:  geometry readout\n\
   Click on a trajectory:  particle & trajectory readout\n\
   Ctrl + click on a volume:  see daughters.\n\
   Shift + click on a volume:  see mother.\n\n\
EXTENDED PICKING mode (Arrow+ viewer button):\n\n\
   Hover the mouse over a volume or trajectory for\n\
   overlayed readout.\n\n\
ELEMENT NAVIGATION (requires Reference Path):\n\n\
   Click on element in list:  centers view on element\n\
   Arrow keys:  rotate in 90 degree steps around element  \n\
   Shift + Right Arrow:  move to next element\n\
   Shift + Left Arrow:  move to previous element\n\n\
FLY mode (requires Reference Path):\n\n\
   Page Up:  Increase speed\n\
   Page Down:  Decrease speed (& reverse if wanted)\n\
   Up Arrow:  raise camera above path\n\
   Down Arror:  lower camera below path\n\
   Escape:  Exit fly mode";
   helpmsgbox->setText(messagetxt);
   helpmsgbox->setModal(false);
   //   helpmsgbox->setWindowModality(Qt::NonModal);

   // AUXILIARY LISTS WINDOW

   // Bypass the namespace in order to make a persistent object
   AuxWindowDialog = new Ui_Dialog;
   AuxWindow = new QDialog(parent);
   AuxWindowDialog->setupUi(AuxWindow);

   // SIGNALS
   connect(AuxWindowDialog->listWidget, SIGNAL(itemClicked(QListWidgetItem*)),
           this, SLOT(LoadBookmarkCB(QListWidgetItem*)));
   connect(AuxWindowDialog->listWidget1, SIGNAL(itemClicked(QListWidgetItem*)),
           this, SLOT(LookAtSceneElementCB(QListWidgetItem*)));
   connect(AuxWindowDialog->pushButton_2, SIGNAL(clicked()),
           this, SLOT(DeleteBookmarkCB()));
   connect(AuxWindowDialog->pushButton_3, SIGNAL(clicked()),
           this, SLOT(RenameBookmarkCB()));
   connect(AuxWindowDialog->pushButton, SIGNAL(clicked()),
           this, SLOT(SortBookmarksCB()));
   
   // FWJ Better to do this after viewer window is realized
   //   AuxWindow->show();
   //   AuxWindow->raise();
   //   AuxWindow->activateWindow();
}


// Called right after buttons and widgets get realized.
// It sets the viewpoint last accessed.
void G4OpenInventorQtExaminerViewer::afterRealizeHook()
{
   SoQtExaminerViewer::afterRealizeHook();

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
         QMessageBox msgbox;
         msgbox.setFont(*font);
         QString messagetxt = "Error reading bookmark file ";
         messagetxt.append(QString(fileName.c_str()));
         msgbox.setText(messagetxt);
         msgbox.exec();
      } else {
         // Opens a file without erasing it
         fileOut.open(fileName.c_str(), std::ios::in);
         fileOut.seekp(0, std::ios::end); // For appending new data to the end
         // FWJ DEBUG
         // G4cout << "afterRealizeHook: opened EXISTING bookmark file"
         //        << G4endl;
         if (viewPtList.size()) {
            // FWJ disabled auto-selection of first viewpoint.
            // Initial view should be user-controllable & not forced
            //    setViewPt();
            addViewPoints();
         }
      }
      fileIn.close();
   } else {
      // Creates a new default bookmark file
      fileOut.open(fileName.c_str());
      // FWJ DEBUG
      // G4cout << "afterRealizeHook: Opened a NEW bookmark file" << G4endl;
   }

   fileIn.clear();

   SoSeparator* root = (SoSeparator*) (getSceneManager()->getSceneGraph());
   if (root == NULL)
      SoDebugError::post("G4OpenInventorQtExaminerViewer::afterRealizeHook", "Root is null.");
   else {
      root->addChild(myCam); // For position/orientation calculation during animation
   }

   sceneChangeSensor = new SoNodeSensor;
   sceneChangeSensor->setFunction(sceneChangeCB);
   sceneChangeSensor->attach(root);
   sceneChangeSensor->setData(this);

   ///////////////////////////// MOUSEOVER & PICK /////////////////////

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

   ///////////////////////////// MOUSEOVER & PICK /////////////////////

   AuxWindow->show();
   AuxWindow->raise();
   AuxWindow->activateWindow();

   auto UI = G4UImanager::GetUIpointer();
   uiQt = dynamic_cast<G4UIQt*>(UI->GetG4UIWindow());
   // This explicitly sets the TabWidget as parent before addTab():
   if (uiQt) {
      viewerParent = getParentWidget();
      viewerParent2 = viewerParent->parentWidget();
      uiQt->AddTabWidget(getParentWidget(), *fName);
      uiQtTabIndex = uiQt->GetViewerTabWidget()->currentIndex();
      //      attached = TRUE;
      addAppPushButton(detachButton);
   }
}


// This method locates a named node in the superimposed or original scene.
// FWJ RENAME THIS
SoNode*
G4OpenInventorQtExaminerViewer::getSuperimpositionNode(SoNode* root,
                                                       const char* name)
{
   if (!searcher)
      searcher = new SoSearchAction;
   searcher->reset();
   searcher->setName(SbName(name));
   searcher->setInterest(SoSearchAction::FIRST);
   searcher->setSearchingAll(TRUE);
   searcher->apply(root);
   assert(searcher->getPath());
   return searcher->getPath()->getTail();
}


// FWJ don't know why userdata is called "closure"
// It contains the this pointer!
void G4OpenInventorQtExaminerViewer::superimpositionCB(void * closure,
                                                       SoAction * action)
{
   if (closure)
      ((G4OpenInventorQtExaminerViewer*)closure)->superimpositionEvent(action);
}


// Renders and positions speed indicator and longitudinal
// distance/viewpoint name on the drawing canvas
void G4OpenInventorQtExaminerViewer::superimpositionEvent(SoAction * action)
{

   if (!action->isOfType(SoGLRenderAction::getClassTypeId()))
      return;
   SbViewportRegion vpRegion =
      ((SoGLRenderAction*)action)->getViewportRegion();
   SbVec2s viewportSize = vpRegion.getViewportSizePixels();
   
   // Aspect is WIDTH/HEIGHT
   float aspect = float(viewportSize[0]) / float(viewportSize[1]);

   // FWJ DEBUG
   //   G4cout << "SPEVENT X0 Y0 DX DY aspect: " << vpRegion.getViewportOrigin()[0] <<
   //      " " << vpRegion.getViewportOrigin()[1] <<
   //      " " << viewportSize[0] <<
   //      " " << viewportSize()[1] <<
   //      " " << aspect << G4endl;

   // Translation and scale factor for animation speed indicator...

   float factorx = 1.0f / float(viewportSize[1]) * 220.0f;
   float factory = factorx;

   if (aspect > 1.0f) {
      stranslation->translation.setValue(SbVec3f(0.0f, -0.4f, 0.0f));
   } else {
      stranslation->translation.setValue(SbVec3f(0.0f, -0.4f / aspect, 0.0f));
      factorx /= aspect;
      factory /= aspect;
   }
   if (viewportSize[0] > 500)
      factorx *= 500.0f / 400.0f;
   else
      factorx *= float(viewportSize[0]) / 400.0f;

   sscale->scaleFactor.setValue(SbVec3f(factorx, factory, 1.0f));

   // TEXT OVERLAY...

   // FWJ Simplified and rewrote the following section to ease problems
   // with the overlayed text after a viewer window resize.
   // Result is now readable but needs further refinement of the scaling.

   float xInfo, yInfo, xLogName, yLogName, xSolid, ySolid,
      xMaterial, yMaterial, xZPos, yZPos;

   // Base point for navigation distance or viewpoint name
   // Origin is at center of render area.
   xInfo = -.475;
   yInfo = .475;
   // Menu bar height in same coordinates:
   float mbgap = 0.03;
   if (aspect > 1.) xInfo = xInfo*aspect;
   if (aspect < 1.) yInfo = yInfo/aspect;
   yInfo = yInfo - mbgap*aspect;

   // Following are relative to above base point
   xLogName = 0.0;
   yLogName = -.88 + mbgap*aspect;
   xSolid = 0.0;
   ySolid = -.91 + mbgap*aspect;
   xMaterial = 0.0;
   yMaterial = -.94 + mbgap*aspect;
   xZPos = 0.0;
   yZPos = -.97 + mbgap*aspect;

   // Top line
   curInfoTrans->translation.setValue(SbVec3f(xInfo, yInfo, 0.0));

   // Bottom lines
   mouseOverTransLogName->translation.setValue(SbVec3f(xLogName, yLogName, 0.0));
   mouseOverTransSolid->translation.setValue(SbVec3f(xSolid, ySolid, 0.0));
   mouseOverTransMaterial->translation.setValue(SbVec3f(xMaterial, yMaterial, 0.0));
   mouseOverTransZPos->translation.setValue(SbVec3f(xZPos, yZPos, 0.0));

   if (currentState == VIEWPOINT) { // Displaying viewpoint name
      curInfoFont->size.setValue(15);
      curInfoFont->name.setValue("defaultFont:Italic");
      curInfoText->string.setValue(SbString(curViewPtName));
   }
   else if(currentState == GENERAL) { // Displaying longitudinal distance
      curInfoFont->size.setValue(16);
      curInfoFont->name.setValue("defaultFont:Bold");
      curInfoText->string.setValue(SbString(""));
   }
   else {
      if (refParticleIdx < (int) refParticleTrajectory.size() - 1) {
         curInfoFont->size.setValue(16);
         curInfoFont->name.setValue("defaultFont:Bold");
         char zPos[20];
         // FWJ need a better format here
         snprintf(zPos, sizeof zPos, "%-7.2f [m]", refZPositions[refParticleIdx] / 1000);
         curInfoText->string.setValue(SbString(zPos));
      }
   }
}


//  Loads view point data from a file into a vector.

bool G4OpenInventorQtExaminerViewer::loadViewPts() 
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

      std::size_t end = token.find_last_not_of(' '); // Remove padded spaces
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


// Rotates camera 90 degrees around a scene element.
// Rotation is animated for smoothness.
void G4OpenInventorQtExaminerViewer::rotateCamera()
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
void G4OpenInventorQtExaminerViewer::moveCamera(float dist, bool lookdown)
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
   else {

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
         refParticleIdx = (int) refParticleTrajectory.size() - 2;
      }

      // Set start and end points
      p1 = refParticleTrajectory[refParticleIdx];
      p2 = refParticleTrajectory[refParticleIdx + step];

      // Get the direction from p1 to p2
      particleDir = p2 - p1;
      particleDir.normalize();

      if(prevParticleDir == SbVec3f(0,0,0)) {
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


void G4OpenInventorQtExaminerViewer::pickingCB(void *aThis, 
                                               SoEventCallback *eventCB)
{
   SoHandleEventAction* action = eventCB->getAction();
   const SoPickedPoint *pp = action->getPickedPoint();
   G4OpenInventorQtExaminerViewer* This = (G4OpenInventorQtExaminerViewer*)aThis;

   if(pp != NULL) {

      SoPath* path = pp->getPath();
      SoNode* node = ((SoFullPath*)path)->getTail();

      if(node->getTypeId() == SoLineSet::getClassTypeId()) {

         if(This->pickRefPathFlag) {
            This->pickRefPathFlag = false;
            if(This->viewingBeforePickRef != This->isViewing())
               This->setViewing(This->viewingBeforePickRef);
            else
               This->setComponentCursor(SoQtCursor(SoQtCursor::DEFAULT));

            // The trajectory is a set of lines stored in a LineSet
            SoLineSet * trajectory = (SoLineSet *)node;
            // FWJ DEBUG
            // G4cout << "FOUND trajectory LineSet" << trajectory << G4endl;

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
               if(tmpNode->getTypeId() == SoCoordinate3::getClassTypeId()) {
                  //node found
                  coords = (SoCoordinate3 *)tmpNode;
                  break;
               }
            }

            if(coords == NULL) {
               G4warn << "Could not find the coordinates node"
                  " for the picked trajectory." << G4endl;
               G4warn << " Reference trajectory not set" << G4endl;
               return;
            }
            // FWJ DEBUG
            // G4cout << "FOUND SoCoordinate3 node " << coords << G4endl;


            if ((This->lshiftdown)	|| (This->rshiftdown))
               This->setReferencePath(trajectory, coords, true);  //APPENDING
            else
               This->setReferencePath(trajectory, coords, false);

            return;

         }
         else if(This->abbrOutputFlag) {

            G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
            if(attHolder && attHolder->GetAttDefs().size()) {

               std::string strTrajPoint = "G4TrajectoryPoint:";
               std::ostringstream oss;
               for (std::size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
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
               G4warn << "SoNode : " << node
                      << " SoType : " << cls
                      << " name : " << name
                      << G4endl;
               G4warn << "No attributes attached." << G4endl;
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
         for (std::size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
            G4cout << G4AttCheck(attHolder->GetAttValues()[i],
                                 attHolder->GetAttDefs()[i]);
         }
      } else {
         G4String name((char*)node->getName().getString());
         G4String cls((char*)node->getTypeId().getName().getString());
         G4warn << "SoNode : " << node
                << " SoType : " << cls
                << " name : " << name
                << G4endl;
         G4warn << "No attributes attached." << G4endl;
      }

      //Suppress other event handlers
      eventCB->setHandled();
   }
}


void G4OpenInventorQtExaminerViewer::mouseoverCB(void *aThis, SoEventCallback *eventCB)
{
   SoHandleEventAction* action = eventCB->getAction();
   const SoPickedPoint* pp = action->getPickedPoint();
   G4OpenInventorQtExaminerViewer* This = (G4OpenInventorQtExaminerViewer*)aThis;

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
         SbVec3f centr = bBox.getCenter();
         centr.getValue(x,y,z);
         ssZPos << "Pos:  " << x << "  " << y << "  " << z;

         G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
         if(attHolder && attHolder->GetAttDefs().size()) {

            std::vector<const std::map<G4String,G4AttDef>*> vecDefs =
               attHolder->GetAttDefs();
            std::vector<const std::vector<G4AttValue>*> vecVals =
               attHolder->GetAttValues();
            for (std::size_t i = 0; i < vecDefs.size(); ++i) {
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
         // G4cout << "Trajectory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
         G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(node);
         if(attHolder && attHolder->GetAttDefs().size()) {
            std::string strTrajPoint = "G4TrajectoryPoint:";
            std::ostringstream oss;
            G4String t1, t1Ch, t2, t3, t4;
            for (std::size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
               // G4cout << "Getting index " << i << " from attHolder" << G4endl;
               // No, returns a vector!
               //   G4AttValue* attValue = attHolder->GetAttValues()[i];
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
                     std::size_t ipos = value.rfind(" ");
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
//             G4cout << "  NOW CALLING G4AttCheck" << G4endl;
//             G4cout << G4AttCheck(attHolder->GetAttValues()[i],
//                                     attHolder->GetAttDefs()[i]);
//             oss << G4AttCheck(attHolder->GetAttValues()[i],
//                                  attHolder->GetAttDefs()[i]);
//             if(oss.str().find(strTrajPoint) != std::string::npos) {
//                // Last attribute displayed was a trajectory point.  Since we
//                // want abbreviated output, display the last one and exit
//                // (unless we're already at the last (and only) trajectory point)
//                if(i != attHolder->GetAttDefs().size()-1) {
//                   G4cout << G4AttCheck(
//                      attHolder->GetAttValues()[attHolder->GetAttDefs().size()-1],
//                      attHolder->GetAttDefs()[attHolder->GetAttDefs().size()-1]);
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
      if(std::string(This->mouseOverTextMaterial->string.getValues(0)->getString()) != ssMaterials.str()) {
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


// Called by hitting PageUp during animation.
void G4OpenInventorQtExaminerViewer::incSpeed() {
   if (std::ceil(animateBtwPtsPeriod * 100) >= 4) {
      if (speedStep > 0.08)
         speedStep -= 0.02;
      else
         speedStep = 0.02;
      animateBtwPtsPeriod -= speedStep;
   } else
      animateBtwPtsPeriod = 0.0;

   if (currentState != PAUSED_ANIMATION) {
      int lastIdx = (int) refParticleTrajectory.size() - 1;
      if (refParticleIdx < lastIdx && !animateSensor->isScheduled())
         animateRefParticle();
   }
}

// Called by hitting PageDown during animation.
void G4OpenInventorQtExaminerViewer::decSpeed() {
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


// Based on the user's interaction the speed indicator bar needs to be adjusted

void G4OpenInventorQtExaminerViewer::updateSpeedIndicator(void)
{
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


void G4OpenInventorQtExaminerViewer::actualRedraw(void) {
	switch (currentState) {
	case ANIMATION:
	case REVERSED_ANIMATION:
	case PAUSED_ANIMATION:
		updateSpeedIndicator();
		SoQtExaminerViewer::actualRedraw();
		break;
	default:
		SoQtExaminerViewer::actualRedraw();
		break;
	}
}


void G4OpenInventorQtExaminerViewer::setReferencePath(SoLineSet *lineset,
       SoCoordinate3 *coords, bool append)
{
   // TODO:  Color the reference path
   // Disable the color stuff for now: changes all trajectories
   // FWJ See G4OpenInventorXtExaminerViewer.cc for test code

   // The trajectory is composed of all the polyline segments in the
   // multiple value field (SoMFInt32) numVertices.
   // For each of the numVertices.getNum()* polyline segments,
   // retrieve the points from the SoCoordinate3 node

   SbVec3f refParticlePt;

   if(!append)
      refParticleTrajectory.clear();

   for(int i = 0; i < lineset->numVertices.getNum(); ++i) {
      for(int j = 0; j < lineset->numVertices[i]; ++j) {
         refParticlePt = coords->point[j];
         refParticleTrajectory.push_back(refParticlePt);
      }
   }
   // Remove points that are too close to each other
   evenOutRefParticlePts();
   setReferencePathZPos();
   getSceneElements();
   sortElements();
}


void G4OpenInventorQtExaminerViewer::setReferencePathZPos()
{
   refZPositions.clear();
   refZPositions.push_back(0);
   float dist;
   for(unsigned int i=0; i < refParticleTrajectory.size() - 1; ++i) {
      dist = (refParticleTrajectory[i] - 
              refParticleTrajectory[i + 1]).length();
      refZPositions.push_back(refZPositions[i] + dist);
   }
}


void G4OpenInventorQtExaminerViewer::findAndSetRefPath()
{
   SoSearchAction action;
   action.setType(SoLineSet::getClassTypeId(),false);
   action.setInterest(SoSearchAction::ALL);
   action.apply(getSceneGraph());

   SoPathList &pathList = action.getPaths();

   if(pathList.getLength() != 0) {

      SoCoordinate3 * coords = NULL;
      std::vector<SoCoordinate3 *> coordvec;
      std::vector<SoLineSet *> linevec;

      bool refPathFound = false;
      for(int i = 0; i < pathList.getLength(); ++i) {
         SoFullPath *path = (SoFullPath *)pathList[i];

         G4AttHolder* attHolder = dynamic_cast<G4AttHolder*>(path->getTail());
         for (std::size_t j = 0; j < attHolder->GetAttDefs().size(); ++j) {
            std::ostringstream oss;
            oss << G4AttCheck(attHolder->GetAttValues()[j],
                              attHolder->GetAttDefs()[j]);

            std::string findStr = "Type of trajectory (Type): ";
            std::string compareValue = "REFERENCE";
            std::size_t idx = oss.str().find(findStr);

            if(idx != std::string::npos) {
               if(oss.str().substr(idx + findStr.size(),
                                   compareValue.size()) == compareValue) {
                  coords = getCoordsNode(path);
                  if(coords != NULL) {
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
                     break;	//Not a primary track, continue with next track

                  coords = getCoordsNode(path);
                  if(coords != NULL) {
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

      if(refPathFound) {
         //set ref path to last traj, coord in the vecs
         setReferencePath(linevec.back(), coordvec.back());
         return;
      }
      //else

      int longestIdx = 0;
      float longestLength = 0.0;
      // For all paths
      for(unsigned int i=0;i < linevec.size(); ++i) {

         //First generate a vector with all the points in this lineset
         std::vector<SbVec3f> trajectory;
         // For all lines in the i path
         for(int j=0; j < linevec[i]->numVertices.getNum(); ++j) {
            // For all points in line j
            for(int k=0; k < linevec[i]->numVertices[j]; ++k) {
               trajectory.push_back(coordvec[i]->point[k]);
            }
         }

         // Then calculate the total length
         float tmpLength=0.0;
         for(unsigned int j=0; j < trajectory.size() - 1; ++j) {
            tmpLength += (trajectory[j] - trajectory[j + 1]).length();
         }

         if(tmpLength > longestLength) {
            longestIdx = i;
            longestLength = tmpLength;
         }
      }

      // Set the longest path as the reference path
      setReferencePath(linevec[longestIdx], coordvec[longestIdx]);
   }
}


SoCoordinate3 * G4OpenInventorQtExaminerViewer::getCoordsNode(SoFullPath *path)
{
   SoLineSet *trajectory = (SoLineSet *)path->getTail();
   SoSeparator * grpNode = (SoSeparator*)(((SoFullPath*)path)->getNodeFromTail(1));
   int nodeIndex = grpNode->findChild(trajectory);
   SoNode * tmpNode;

   // We allow only 100 iterations, in case the node isn't found
   // (should take only a few iterations)
   for (int i = 0; i < 100; ++i) {
      --nodeIndex;

      tmpNode = grpNode->getChild(nodeIndex);
      if(tmpNode->getTypeId() == SoCoordinate3::getClassTypeId()) {
         //node found
         return (SoCoordinate3 *)tmpNode;
      }
   }
   return NULL;	//coords node not found
}


// Displays scene elements on the right side of listsDialog.
// else: scene graph is searched for Geant4_SoPolyhedron type nodes
void G4OpenInventorQtExaminerViewer::getSceneElements()
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

   // FWJ DEBUG
   //   G4cout << "Searching for elements....." << G4endl;
   search.apply(root);

   SoPathList &pl = search.getPaths();


   // First find which names occur more than once so we can append a counter to them
   for (int i = 0; i < pl.getLength(); i++) {
      SoFullPath *path = (SoFullPath *)pl[i];
      node = (Geant4_SoPolyhedron *)path->getTail();
      eltName = node->getName();
      //      G4cout << "  FOUND " << i << "  " << eltName << G4endl;
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
      sceneElements.push_back(el);
   }
}


float G4OpenInventorQtExaminerViewer::sqrlen(const SbVec3f &a)
{
   float x,y,z;
   a.getValue(x,y,z);
   return x*x + y*y + z*z;
}


void G4OpenInventorQtExaminerViewer::distanceToTrajectory(const SbVec3f &q,
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
   // 				p(t) = a+t?(b?a)
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

   const std::size_t count = refParticleTrajectory.size();
   assert(count>0);

   SbVec3f b = refParticleTrajectory[0];
   SbVec3f dbq = b - q;
   float sqrDist = sqrlen(dbq);
   closestPoint = b;
   index = 0;
   for (std::size_t i = 1; i < count; ++i) {
      const SbVec3f a = b;
      const SbVec3f daq = dbq;
      b = refParticleTrajectory[i];
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

      if (t<0.) {
         // The trajectory point occurs before point a
         // Go to the next point
         continue;
      }
      float current_dist;
      if (t<=1.) {
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

      if (current_dist < sqrDist) {
         sqrDist = current_dist;
         closestPoint = a + t*(b-a);
         index = (int) i;
      }
   }

   dist = std::sqrt(sqrDist);
}


void G4OpenInventorQtExaminerViewer::sortElements()
{
   if(refParticleTrajectory.empty())
      return;

   float * trajLength = new float[refParticleTrajectory.size()];
   typedef std::map<elementForSorting, sceneElement> sortedMap;
   sortedMap sorted;

   // For every point on the reference trajectory, compute
   // the total length from the start
   SbVec3f prevPoint;
   std::vector<SbVec3f>::iterator itRef = refParticleTrajectory.begin();
   int trajIndex = 0;
   prevPoint = *itRef;
   trajLength[trajIndex] = 0.0;
   ++itRef;
   ++trajIndex;
   for(; itRef != refParticleTrajectory.end(); ++itRef, ++trajIndex) {
      trajLength[trajIndex] = trajLength[trajIndex-1] + (*itRef - prevPoint).length();
      prevPoint = *itRef;
   }

   // Compute the smallest distance between the element
   // and the reference trajectory (find the closest point),
   // then map the element to the trajectory length of that
   // point (calculated above)
   SoGetBoundingBoxAction bAction(getViewportRegion());
   SbVec3f elementCoord;
   std::vector<sceneElement>::iterator itEl;
   int elementIndex;
   elementForSorting el;
   for(itEl = sceneElements.begin(), elementIndex = 0;
       itEl != sceneElements.end(); ++itEl, ++elementIndex) {
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
      el.distanceToBeamlineStart = (itEl->center - refParticleTrajectory[0]).length();

      // This map of the scene elements (or their coordinates rather)
      // is automatically sorted by trajectory length (Z coord), then
      // by the distance between the element and the point in case the Z coord
      // is the same as another element.  This is done by using as a key
      // an element structure which implements the operator for weak ordering
      sorted.insert(std::make_pair(el,*itEl));
   }

   // store the sorted elements into the vector field
   sceneElements.clear();

   sortedMap::iterator itSorted = sorted.begin();
   for(; itSorted != sorted.end(); itSorted++)
      sceneElements.push_back(itSorted->second);

   zcoordSetFlag = true;

   createElementsList();

   delete[] trajLength;
}


void G4OpenInventorQtExaminerViewer::createElementsList()
{
   // FWJ DEBUG
   //   G4cout << "Populating ELEMENT LIST..." << G4endl;

   AuxWindowDialog->listWidget1->clear();
   //   int size = sceneElements.size();

   std::vector<sceneElement>::const_iterator it;
   std::stringstream ss;

   for(it=sceneElements.begin(); it!=sceneElements.end(); ++it) {
      ss << it->name;
      if(zcoordSetFlag)
         ss << " [" << it->closestPointZCoord << "]";

      new QListWidgetItem(ss.str().c_str(), AuxWindowDialog->listWidget1); 
      ss.str("");
   }
}


// Called when user clicks a scene element in listsDialog.
// Zooms onto that element.
void
G4OpenInventorQtExaminerViewer::LookAtSceneElementCB(QListWidgetItem* item)
{
   char* value;
   std::string elementField;

   // FWJ DEBUG
   //   G4cout << "AuxWindow: listWidget1 select element CALLBACK" << G4endl;

   SoCamera * cam = getCamera();

   if (SoQtExaminerViewer::isAnimating())
      stopAnimating();

   value = strdup(qPrintable(item->text()));
   //   G4cout << "LOOKING FOR BOOKMARK " << value << G4endl;

   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == PAUSED_ANIMATION ) {
      if (animateSensor->isScheduled())
         animateSensor->unschedule();
      setSuperimpositionEnabled(superimposition, FALSE);
      maxSpeed = 0.0f;
      scheduleRedraw();
      restoreCamera();
      currentState = prevState;
   } else if (currentState == VIEWPOINT)
      setSuperimpositionEnabled(superimposition, FALSE);

   elementField = value;

   std::size_t idx = elementField.find_last_of("[");
   if(idx == std::string::npos)
      idx = elementField.size(); //if "[" not found for whatever reason (list not sorted)
   else
      idx--; // To get rid of the space that is between the name and '['

   bool error = false;
   SoFullPath *path;
   SoSearchAction search;
   SoNode *root = getSceneManager()->getSceneGraph();
   int counter;
   std::size_t idxUnderscore = elementField.find_last_of("_");

   parseString<int>(counter, 
                          elementField.substr(idxUnderscore + 1, idx), error);

   SoBaseKit::setSearchingChildren(TRUE);
   search.reset();
   search.setSearchingAll(TRUE);

   // G4cout << "  Starting search for elementField " << elementField 
   //        << G4endl;

   if(error) { // No counter is present => element name was not modified
      curEltName = elementField.substr(0, idx);
      search.setName(curEltName.c_str());
      search.apply(root);

      path = (SoFullPath *)search.getPath();
   }
   else {
      curEltName = elementField.substr(0, idxUnderscore);
      search.setInterest(SoSearchAction::ALL);
      search.setName(curEltName.c_str());
      search.apply(root);

      SoPathList &pl = search.getPaths();
      path = (SoFullPath *)pl[counter - 1]; // Since counter starts at 1, not 0
   }

   G4ThreeVector global;

   // FWJ FLIP THIS
   if ((idx > 0) && (path)) {

      if(!refParticleTrajectory.empty()) {

         SoGetBoundingBoxAction bAction(getViewportRegion());
         bAction.apply(path);
         SbBox3f bBox = bAction.getBoundingBox();
         SbVec3f elementCoord = bBox.getCenter();

         refParticleIdx = 0;
         SbVec3f p;

         float absLengthNow, absLengthMin;
         int maxIdx = (int) refParticleTrajectory.size() - 2;
         int targetIdx = 0;
         SbVec3f dir;

         p = refParticleTrajectory[refParticleIdx];
         absLengthMin = (p - elementCoord).length();
         refParticleIdx++;

         // Find a ref. particle's point closest to element's global coords
         while (refParticleIdx < maxIdx) {
            p = refParticleTrajectory[refParticleIdx];
            absLengthNow = (p - elementCoord).length();

            if (absLengthNow < absLengthMin) {
               absLengthMin = absLengthNow;
               targetIdx = refParticleIdx;
            }
            refParticleIdx++;
         }

         if (currentState != BEAMLINE) { // Set up default zoom
            SbVec3f p1, pN;
            currentState = BEAMLINE;
            prevParticleDir = SbVec3f(0,0,0); //so that moveCamera() knows sets default parameters
            
            p1 = prevPt = refParticleTrajectory[0];
            pN = refParticleTrajectory[refParticleTrajectory.size() - 1];
            distance = (pN - p1).length() / 10;

            // FWJ Rather than switching to a default height, it is more flexible
            // to keep the same height(magnification) while moving the camera.
            // if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
            //    ((SoOrthographicCamera *) cam)->height.setValue(defaultHeight);
            // // FWJ Restore the default height instead of hard-wired value
            // // ((SoOrthographicCamera *) cam)->height.setValue(10000.0f);
            // }
            // else if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))

            // FWJ required to avoid extreme perspective after camera move:
            if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))
               ((SoPerspectiveCamera*)cam)->heightAngle.setValue(defaultHeightAngle);

         } else {
            if (cam->isOfType(SoPerspectiveCamera::getClassTypeId()))
               distance = (prevPt - cam->position.getValue()).length();
         }
         refParticleIdx = targetIdx;

         //////////////////////////////////////////////////////////////
         setSuperimpositionEnabled(superimposition, TRUE);
         axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
         animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
         animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);
         scheduleRedraw();
         //////////////////////////////////////////////////////////////

         moveCamera(distance);

      }
      
      else {
         offsetFromCenter.setValue(0, 0, 1);
         distance = 50;// small number since using viewAll() for default zoom
         upVector.setValue(0, 1, 0);
         moveCamera(distance);
         cam->viewAll(path, getViewportRegion());
      }
   }

}


void G4OpenInventorQtExaminerViewer::FileLoadRefPathCB()
{
   //   G4cout << "File: Load Ref Path CALLBACK" << G4endl;

   QFileDialog filedialog(getParentWidget(), tr("Load Reference Path"));
   filedialog.setFileMode(QFileDialog::AnyFile);
   filedialog.setFont(*font);
   if (!filedialog.exec()) return;
   QStringList filenameinlist = filedialog.selectedFiles();
   QString filenamein = filenameinlist[0];

   //   G4cout << "Input file name is " << qPrintable(filenamein) << G4endl;

   char* filename = new char[filenamein.size()+1];
   filename = strdup(qPrintable(filenamein));
   //   G4cout << "char[] file name is " << filename << G4endl;

   std::ifstream ifs(filename);
   if(ifs.is_open()) {
      refParticleTrajectory.clear();
      float x,y,z;
      while(ifs >> x >> y >> z) {
         refParticleTrajectory.push_back(SbVec3f(x,y,z));
      }
      ifs.close();
   } else {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Reference Path file not found: ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      return;
   }
   if (refParticleTrajectory.size() < 2) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Invalid Reference Path";
      msgbox.setText(messagetxt);
      msgbox.exec();
      return;
   }
   // Following setReferencePath() ...
   evenOutRefParticlePts();
   setReferencePathZPos();
   getSceneElements();
   sortElements();
}


void G4OpenInventorQtExaminerViewer::FileSaveRefPathCB()
{
   //   G4cout << "File: Save Ref Path CALLBACK" << G4endl;

   QFileDialog filedialog(getParentWidget(), tr("Save Reference Path"));
   filedialog.setFileMode(QFileDialog::AnyFile);
   // To enable confirmation of overwriting
   filedialog.setAcceptMode(QFileDialog::AcceptSave);
   filedialog.setFont(*font);
   if (!filedialog.exec()) return;
   QStringList filenameinlist = filedialog.selectedFiles();
   QString filenamein = filenameinlist[0];

   //   G4cout << "Input file name is " << qPrintable(filenamein) << G4endl;

   char* filename = new char[filenamein.size()+1];
   filename = strdup(qPrintable(filenamein));
   //   G4cout << "char[] file name is " << filename << G4endl;

   std::ofstream ofs(filename);
   if (ofs.is_open()) {
      float x,y,z;
      for (unsigned int i=0; i < refParticleTrajectory.size(); ++i) {
         refParticleTrajectory[i].getValue(x,y,z);
         ofs << x << " " << y << " " << z << "\n";
      }
      ofs.close();
   } else {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error opening file ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
   }

}

void G4OpenInventorQtExaminerViewer::evenOutRefParticlePts()
{
   if(refParticleTrajectory.empty())
      return;

   SbVec3f p1, p2, p3, dirNow, dirNxt, dir, p2_tmp, p_start, p_corner, p_nxt;
   float avgDistBtwPts = 0;
   float totalDistBtwPts = 0;
   std::vector<SbVec3f> newRefParticleTrajectory;
   SbVec3f refPoint;
   std::size_t size = refParticleTrajectory.size() - 1;
   int numOfPts = 0;
   for (std::size_t i = 0; i < size; ++i) {
      p1 = refParticleTrajectory[i];
      p2 = refParticleTrajectory[i + 1];
      if (p1 == p2)
         continue;
      numOfPts++;
      totalDistBtwPts += (p2 - p1).length();
   }
   // Nothing useful to do (and fix Coverity)
   if (numOfPts <= 2) return;

   avgDistBtwPts = totalDistBtwPts / numOfPts;
   float minDistAllowed = 0.75 * avgDistBtwPts;
   //	float maxDistAllowed = 1.25 * avgDistBtwPts; // Pts tend to be close not far

   float x, y, z;
   std::size_t i = 0, j = 0;
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


void G4OpenInventorQtExaminerViewer::saveCurCamera()
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


void G4OpenInventorQtExaminerViewer::restoreCamera()
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


void G4OpenInventorQtExaminerViewer::animateSensorRotationCB(void *data, 
                                                             SoSensor *sensor)
{
   SbTime curTime = SbTime::getTimeOfDay();
   G4OpenInventorQtExaminerViewer* This = (G4OpenInventorQtExaminerViewer*) data;

   SoTimerSensor* s = (SoTimerSensor*) sensor;

   float t = float((curTime - s->getBaseTime()).getValue())
      / This->animateBtwPtsPeriod;

   if ((t > 1.0f) || (t + s->getInterval().getValue() > 1.0f))
      t = 1.0f;
   SbBool end = (t == 1.0f);

   if (end) {
      This->animateSensorRotation->unschedule();
      if(This->rotCnt) {
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

void G4OpenInventorQtExaminerViewer::animateSensorCB(void *data, 
                                                     SoSensor *sensor)
{
   SbTime curTime = SbTime::getTimeOfDay();
   G4OpenInventorQtExaminerViewer* This = (G4OpenInventorQtExaminerViewer*) data;
   SoCamera *cam = This->getCamera();
   SoTimerSensor* s = (SoTimerSensor*) sensor;

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


void G4OpenInventorQtExaminerViewer::setStartingPtForAnimation()
{
   if (SoQtExaminerViewer::isAnimating())
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


void G4OpenInventorQtExaminerViewer::gotoRefPathStart()
{
   G4OpenInventorQtExaminerViewer::ToolsRefPathStartCB();
}


void G4OpenInventorQtExaminerViewer::ToolsRefPathStartCB()
{
   if (!refParticleTrajectory.size()) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "No current reference path";
      msgbox.setText(messagetxt);
      msgbox.exec();
      return;
   }

   if (currentState == ROTATING)
      return;
   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == PAUSED_ANIMATION) {
      if (animateSensor->isScheduled())
         animateSensor->unschedule();
      setSuperimpositionEnabled(superimposition, FALSE);
      maxSpeed = 0.0f;
      scheduleRedraw();
   } else {
      saveCurCamera();
      prevState = currentState;
      prevRefIdx = refParticleIdx;
   }

   if (SoQtExaminerViewer::isAnimating())
      stopAnimating();

   up_down = 0;
   left_right = 0;
   step = 1;

   refParticleIdx = 0;
   currentState = BEAMLINE;
   setSuperimpositionEnabled(superimposition, TRUE);
   axisSwitch->whichChild.setValue(SO_SWITCH_NONE);
   animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_NONE);
   animSpeedSwitch->whichChild.setValue(SO_SWITCH_NONE);
   scheduleRedraw();

   // FWJ Disabled: this is set in moveCamera()
   // Zoom in
   //   SoCamera *cam = getCamera();
   //   cam->focalDistance = 0.1f;

   prevParticleDir = SbVec3f(0,0,0);

   //Default zoom
   SbVec3f p1 = refParticleTrajectory[0];
   SbVec3f pN = refParticleTrajectory[refParticleTrajectory.size() - 1];
   distance = (pN - p1).length() / 10;

   moveCamera(distance, true);
}


void G4OpenInventorQtExaminerViewer::ToolsRefPathInvertCB()
{
   invertRefPath();
}


void G4OpenInventorQtExaminerViewer::invertRefPath()
{
   std::reverse(refParticleTrajectory.begin(),
                refParticleTrajectory.end());
   setReferencePathZPos();
   sortElements();
}


void G4OpenInventorQtExaminerViewer::animateRefParticle()
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


void G4OpenInventorQtExaminerViewer::addEscapeCallback(void (*callback)())
{
   escapeCallback = callback;
}


void G4OpenInventorQtExaminerViewer::sceneChangeCB(void* userData, SoSensor*)
{
   // FWJ DEBUG
   //   G4cout << "SCENE CHANGE callback" << G4endl;
   // NOTE: could/should be disabled during animation

   G4OpenInventorQtExaminerViewer* This =
      (G4OpenInventorQtExaminerViewer*)userData;
   if(This->newEvents) {
      This->findAndSetRefPath();
      This->newEvents = false;
   }
}


//////////////////////////////////// BOOKMARKS ///////////////////////////

// Adds bookmarks to listsDialog.

void G4OpenInventorQtExaminerViewer::addViewPoints()
{
   std::size_t size = viewPtList.size();
   if (!size) return;

   for (std::size_t i = 0; i < size; ++i) {
      new QListWidgetItem(viewPtList[i].viewPtName,
                          AuxWindowDialog->listWidget); 
   }
}


// Converts a string type word into a float type.

template<class T> 
void G4OpenInventorQtExaminerViewer::parseString(T &t, const std::string &s,
                                                 bool &error) 
{
   std::istringstream str(s);
   if ((str >> t).fail())
      error = true;
}


void
G4OpenInventorQtExaminerViewer::FileOpenBookmarkCB()
{
   // FWJ DEBUG
   //   G4cout << "File: Open Bookmark File CALLBACK" << G4endl;
   QFileDialog filedialog(getParentWidget(), tr("Open bookmark file"));
   filedialog.setFileMode(QFileDialog::ExistingFile);
   filedialog.setFont(*font);
   if (!filedialog.exec()) return;
   QStringList filenameinlist = filedialog.selectedFiles();
   QString filenamein = filenameinlist[0];

   char* filename = new char[filenamein.size()+1];
   filename = strdup(qPrintable(filenamein));
   //   G4cout << "char[] file name is " << filename << G4endl;

   fileIn.close();
   fileIn.open(filename);
   if (fileIn.fail()) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error opening file: ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      //      G4cout << "ERROR opening file " << filename << G4endl;
      fileIn.clear();
      return;
   }
   // Opens a file without erasing it
   cleanUpAfterPrevFile();

   if (!loadViewPts()) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error reading bookmark file: ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      //      G4cout << "ERROR reading bookmark file " << filename << G4endl;
      fileIn.clear();
      return;
   }

   fileName = filename;
   fileOut.open(fileName.c_str(), std::ios::in);
   fileOut.seekp(0, std::ios::end);

   addViewPoints();

   // LATER: display filename in lists window

   fileIn.close();
   fileIn.clear();
}

// Called before loading a new viewpoint file. 
// Resets member fields to default values.

void G4OpenInventorQtExaminerViewer::cleanUpAfterPrevFile()
{
   viewPtIdx = -1;
   viewPtList.clear();
   //   setSuperimpositionEnabled(superimposition, FALSE);
   //   scheduleRedraw();
   currentState = GENERAL;
   if (fileOut.is_open()) fileOut.close();

   AuxWindowDialog->listWidget->clear();   
   AuxWindowDialog->lineEdit->setText(QString(""));
}


void
G4OpenInventorQtExaminerViewer::FileNewBookmarkCB()
{
   //   G4cout << "File: Open New Bookmark File CALLBACK" << G4endl;
   QFileDialog filedialog(getParentWidget(), tr("Open new bookmark file"));
   filedialog.setFileMode(QFileDialog::AnyFile);
   // To enable confirmation of overwriting
   filedialog.setAcceptMode(QFileDialog::AcceptSave);
   // But change the "Save" button text
   filedialog.setLabelText(QFileDialog::Accept, QString("New"));
   filedialog.setFont(*font);
   if (!filedialog.exec()) return;
   QStringList filenameinlist = filedialog.selectedFiles();
   QString filenamein = filenameinlist[0];

   //   G4cout << "Input file name is " << qPrintable(filenamein) << G4endl;

   char* filename = new char[filenamein.size()+1];
   filename = strdup(qPrintable(filenamein));
   //   G4cout << "char[] file name is " << filename << G4endl;

   cleanUpAfterPrevFile();
   fileName = filename;
   fileOut.open(fileName.c_str());
   if (fileOut.fail()) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error opening new bookmark file: ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      // G4cout << "ERROR opening new bookmark file " << filename << G4endl;
   }
}


void
G4OpenInventorQtExaminerViewer::ToolsAnimateRefParticleCB()
{
   //   G4cout << "Tools: Animate Ref Particle CALLBACK" << G4endl;
   if (!refParticleTrajectory.size()) {
      returnToAnim = true;
      G4warn << "No Reference Trajectory" << G4endl;
      return;
   }

   ///////////////////////////////////////////////////////////////
   setSuperimpositionEnabled(superimposition, TRUE);
   maxSpeed = SPEED_INDICATOR_STEP;
   axisSwitch->whichChild.setValue(SO_SWITCH_ALL);
   animSpeedOutlineSwitch->whichChild.setValue(SO_SWITCH_ALL);
   animSpeedSwitch->whichChild.setValue(SO_SWITCH_ALL);
   scheduleRedraw();
   ///////////////////////////////////////////////////////////////

   SoCamera *cam = getCamera();
   //   SbVec3f camDirOld, camDirNew, camDirNew_tmp, camUpVec, P0, P1, P1_tmp;

   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == ROTATING)
      return;

   if (currentState != PAUSED_ANIMATION) {

      saveCurCamera();
      prevState = currentState;
      prevRefIdx = refParticleIdx;

      if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
         toggleCameraType();
         cam = getCamera();
      }

      refParticleIdx = 0; // Set the camera to the starting point of the animation
      animateBtwPtsPeriod = MIN_SPEED;
      speedStep = START_STEP;
      left_right = up_down = 0;

      cam->focalDistance = 0.1f;
      ((SoPerspectiveCamera *) cam)->heightAngle = 0.50f;
   }

   currentState = ANIMATION;
   setStartingPtForAnimation();

   cam->position = (myCam)->position.getValue();
   cam->orientation = (myCam)->orientation.getValue();
   animateRefParticle(); // Animate the camera
}


void
G4OpenInventorQtExaminerViewer::SaveViewPtCB()
{
   //   G4cout << "AppButton: Save Viewpoint CALLBACK" << G4endl;
   // First get viewpoint name ...
   // EMULATING getViewPtNameCB ...
   //   bool ok;
   // Note QString() returns an empty string
   
   // NONE OF THE FOLLOWING CHANGES THE FONT: FORGET IT FOR NOW
   QInputDialog* inputdialog = new QInputDialog(getParentWidget());
   inputdialog->setFont(*font);
   inputdialog->setWindowTitle(tr("Enter a name for the bookmark"));
   inputdialog->setLabelText("Bookmark name");
   //   inputdialog->setTextEchoMode(QLineEdit::Normal);
   inputdialog->adjustSize();
   QString namein;
   if (inputdialog->exec() == QDialog::Accepted)
      namein=inputdialog->textValue().trimmed();
   else
      return;
   if (namein.isEmpty()) return;

   // This easier approach failed: unable to set font size
   //   QString namein = 
   //      QInputDialog::getText(getParentWidget(),
   //                            tr("Enter a name for the bookmark"),
   //                            tr("Bookmark name"), QLineEdit::Normal,
   //                            QString(), &ok);

   const int nVPName = MAX_VP_NAME + 1;
   char* name = new char[nVPName];
   //   strncpy(name, strName.c_str(), nVPName);
   namein.truncate(MAX_VP_NAME);

   QByteArray ba = namein.toLocal8Bit();
   name = strdup(ba.constData());
   // name = strdup(qPrintable(namein))

   // FWJ DEBUG
   //   G4cout << "QString is " << qPrintable(namein) << G4endl;
   //   G4cout << "char[] is  " << name << G4endl;

   for (int i = 0; i < (int)viewPtList.size(); i++) {
      if (!strcmp(name, viewPtList[i].viewPtName)) {
         QMessageBox msgbox;
         msgbox.setText("Bookmark name is already in use");
         msgbox.exec();
         return;
      }
   }

   if (viewPtIdx == -1) viewPtIdx = 0;
   saveViewPt(name);

   saveViewPtItem = new QListWidgetItem(namein,
                                        AuxWindowDialog->listWidget);
   AuxWindowDialog->listWidget->setCurrentItem(saveViewPtItem);
   AuxWindowDialog->lineEdit->setText(namein);
}


// Saves current camera parameters to a viewpoint file.

void G4OpenInventorQtExaminerViewer::saveViewPt(char *name)
{
   SbVec3f axis;
   viewPtData tmp;
   float x, y, z, angle;
   SoCamera* camera = getCamera();

   // NOTE: Xt VSN increments this at end of procedure
   //   viewPtIdx++;

   //  FWJ DEBUG
   //   G4cout << "saveViewPt: saving bookmark " << viewPtIdx << " " << name
   //          << G4endl;

   if (viewPtList.size() == 0) {
      writeViewPtIdx();
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
      SoDebugError::post("G4OpenInventorQtExaminerViewer::saveViewPtCB",
                         "Only Perspective and Orthographic cameras are supported.");
      return;
   }

   viewPtList.push_back(tmp);

   // Now save the view point to a .txt file
   // FWJ DEBUG
   // G4cout << "saveViewPt: writing to Bookmark file " << fileName << G4endl;

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

   // FWJ DEBUG
   //   G4cout << "saveViewPt: finished writing to file" << G4endl <<
   //      "  Next viewPtIdx is " << viewPtIdx << G4endl;
}


// Updates the viewPtIdx in a viewpoint file.

void G4OpenInventorQtExaminerViewer::writeViewPtIdx()
{
   std::string idxStr;
   std::stringstream out;

   out << viewPtIdx;
   idxStr = out.str();
   fileOut.seekp(0, std::ios::beg);

   while ((int) idxStr.length() < MAX_VP_IDX) {
      idxStr += " ";
   }

   // FWJ DEBUG
   //   G4cout << "writeViewPtIdx: " << viewPtIdx << G4endl;
   fileOut << idxStr << "\n";
   fileOut.flush();
   fileOut.seekp(0, std::ios::end);
}


// Receives the name of the bookmark clicked and searches for it in viewPtList.

void G4OpenInventorQtExaminerViewer::LoadBookmarkCB(QListWidgetItem* item)
{
   // FWJ DEBUG
   //   G4cout << "AuxWindow: listWidget LoadBookmark CALLBACK" << G4endl;

   const int nVPName = MAX_VP_NAME + 1;
   char* vpName = new char[nVPName];

   vpName = strdup(qPrintable(item->text()));
   //   G4cout << "LOOKING FOR BOOKMARK " << vpName << G4endl;

   for (int i = 0; i < (int)viewPtList.size(); i++) {
      if (!strcmp(viewPtList[i].viewPtName, vpName)) {
         viewPtIdx = i;
         break;
      }
   }
   //   G4cout << "  FOUND viewPtIdx " << viewPtIdx << G4endl;

   writeViewPtIdx();
   setViewPt();
   AuxWindowDialog->lineEdit->setText(item->text());
}


// Sets the viewpoint based on camera data that viewPtIdx is pointing to.

void G4OpenInventorQtExaminerViewer::setViewPt()
{
   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == ROTATING) {
      if (animateSensor->isScheduled()) animateSensor->unschedule();
      setSuperimpositionEnabled(superimposition, FALSE);
      maxSpeed = 0.0f;
      scheduleRedraw();
   }

   SoCamera * camera = getCamera();
   if (camera == NULL) {
      G4warn << "setViewPt: Camera is null. Unable to set the viewpoint." <<
         G4endl;
      //      String dialogName = (char *) "Missing Camera Node";
      //      std::string msg = "Camera is null. Unable to set the viewpoint.";
      //      warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   if (!viewPtList.size()) {
      G4warn << "setViewPt: There are no viewpoints to load." << G4endl;
      //      String dialogName = (char *) "Missing Viewpoints";
      //      std::string msg = "There are no viewpoints to load.";
      //      warningMsgDialog(msg, dialogName, NULL);
      return;
   }

   if (SoQtExaminerViewer::isAnimating()) stopAnimating();

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
      SoDebugError::post("G4OpenInventorQtExaminerViewer::setViewPt",
                         "Only Perspective and Orthographic cameras are supported.");
      return;
   }

}


void G4OpenInventorQtExaminerViewer::NextViewPtCB()
{
   // FWJ DEBUG
   //   G4cout << "App Button: nextViewPt CALLBACK" << G4endl;

   if (!viewPtList.size()) return;
   if (viewPtIdx >= (int)viewPtList.size() - 1)
      viewPtIdx = 0;
   else
      viewPtIdx++;

   writeViewPtIdx();
   setViewPt();
   char* viewptname = viewPtList[viewPtIdx].viewPtName;
   AuxWindowDialog->lineEdit->setText(QString(viewptname));
}

void G4OpenInventorQtExaminerViewer::PrevViewPtCB()
{
   // FWJ DEBUG
   //   G4cout << "App Button: prevViewPt CALLBACK" << G4endl;

   if (!viewPtList.size()) return;
   if (viewPtIdx == 0)
      viewPtIdx = (int) viewPtList.size() - 1;
   else
      viewPtIdx--;

   writeViewPtIdx();
   setViewPt();
   char* viewptname = viewPtList[viewPtIdx].viewPtName;
   AuxWindowDialog->lineEdit->setText(QString(viewptname));
}


void G4OpenInventorQtExaminerViewer::AbbrOutputCB(bool checked)
{
   // FWJ DEBUG
   //   G4cout << "App Button: abbrOutput CALLBACK" << G4endl;

   abbrOutputFlag = checked;
}


void G4OpenInventorQtExaminerViewer::PickRefPathCB()
{
   // FWJ DEBUG
   //   G4cout << "App Button: pickRefPath CALLBACK" << G4endl;

   // Save viewing state and go to picking mode
   viewingBeforePickRef = isViewing();
   if(isViewing())
      setViewing(false);
   setComponentCursor(SoQtCursor(SoQtCursor::CROSSHAIR));
   pickRefPathFlag = true;
}


void G4OpenInventorQtExaminerViewer::SwitchWireFrameCB(bool checked)
{
   // FWJ DEBUG
   //   G4cout << "App Button: switchWireFrame CALLBACK" << G4endl;

   //   if (switchWireFrameButton->isDown()) {
   if (checked) {
      setDrawStyle(SoQtViewer::STILL, SoQtViewer::VIEW_LINE);
      setDrawStyle(SoQtViewer::INTERACTIVE, SoQtViewer::VIEW_LINE);
   } else {
      setDrawStyle(SoQtViewer::STILL, SoQtViewer::VIEW_AS_IS);
      setDrawStyle(SoQtViewer::INTERACTIVE,
                         SoQtViewer::VIEW_SAME_AS_STILL);
   }
}


void G4OpenInventorQtExaminerViewer::SwitchAxesCB(bool checked)
{
   // FWJ DEBUG
   //   G4cout << "App Button: switchAxes CALLBACK" << G4endl;
   setFeedbackVisibility(checked);
   //   if (checked) {
   //      setFeedbackVisibility(TRUE);
   //   } else {
   //      setFeedbackVisibility(FALSE);
   //   }
}


void G4OpenInventorQtExaminerViewer::DetachCB()
{
   //   FWJ DEBUG
   //   G4cout << "App Button: detach CALLBACK" << G4endl;
   uiQt->GetViewerTabWidget()->removeTab(uiQtTabIndex);
   viewerParent->setParent(viewerParent2);
   removeAppPushButton(detachButton);
   show();
}


void G4OpenInventorQtExaminerViewer::DeleteBookmarkCB()
{
   // FWJ DEBUG
   //   G4cout << "Delete Button: DeleteBookmark CALLBACK" << G4endl;

   // Maybe nothing selected yet
   QListWidgetItem* listitem = AuxWindowDialog->listWidget->currentItem();
   if (!listitem) return;
   if (!(listitem->isSelected())) return;

   QString vpnamein = listitem->text();

   const int nVPName = MAX_VP_NAME + 1;
   char* vpName = new char[nVPName];
   vpName = strdup(qPrintable(vpnamein));
   //   G4cout << "DELETING bookmark " << vpName << G4endl;

   deleteViewPt(vpName);
   delete listitem;
}

// Deletes current viewpoint the user is looking at.
// Updates the input file and bookmarks as well.

void G4OpenInventorQtExaminerViewer::deleteViewPt(char *vpName)
{
   std::string line;
   std::size_t end;
   fileIn.open(fileName.c_str());
   std::ofstream out("temporaryFile.txt");

   if (!vpName)
      vpName = viewPtList[viewPtIdx].viewPtName;

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

   std::size_t idx = 0; // Remove viewpoint from the vector
   std::size_t size = viewPtList.size();
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

   // FWJ check return status: error popups needed here
   int istat = remove(fileName.c_str());
   if (istat == -1) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error removing bookmarks file";
      //      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      //      G4cout << "Error removing bookmarks file" << G4endl;
   }
   istat = rename("temporaryFile.txt", fileName.c_str());
   if (istat == -1) {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error renaming bookmarks file";
      //      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      //      G4cout << "Error renaming bookmarks file" << G4endl;
   }
   fileOut.open(fileName.c_str(), std::ios::in);
   fileOut.seekp(0, std::ios::end);

   if (!viewPtList.size()) { // viewPtList is empty
      curViewPtName = (char *) empty.c_str();
      scheduleRedraw();
   } else {
      if (viewPtIdx >= (int) viewPtList.size())
         viewPtIdx--;
      writeViewPtIdx();
      setViewPt();
   }
}


void G4OpenInventorQtExaminerViewer::RenameBookmarkCB()
{
   // FWJ DEBUG
   //   G4cout << "Rename Button: RenameBookmark CALLBACK" << G4endl;
   // Maybe nothing selected yet
   QListWidgetItem* listitem = AuxWindowDialog->listWidget->currentItem();
   if (!listitem) return;
   if (!(listitem->isSelected())) return;

   QString vpnamein = listitem->text();

   const int nVPName = MAX_VP_NAME + 1;
   //   char* vpName = new char[nVPName];
   //   vpName = strdup(qPrintable(vpnamein));
   //   G4cout << "RENAMING bookmark " << vpName << G4endl;

   QInputDialog* inputdialog = new QInputDialog(getParentWidget());
   inputdialog->setFont(*font);
   inputdialog->setWindowTitle(tr("Enter"));
   inputdialog->setLabelText("New bookmark name");
   inputdialog->adjustSize();
   QString newnamein;
   if (inputdialog->exec() == QDialog::Accepted)
      newnamein=inputdialog->textValue().trimmed();
   else
      return;
   if (newnamein.isEmpty()) return;

   char* newname = new char[nVPName];
   newname = strdup(qPrintable(newnamein));

   std::size_t size = viewPtList.size();
   for (std::size_t i = 0; i < size; ++i) {
      if (!strcmp(newname, viewPtList[i].viewPtName)) {
         QMessageBox msgbox;
         msgbox.setFont(*font);
         msgbox.setText("Bookmark name is already in use");
         msgbox.exec();
      }
   }

   //   G4cout << "RENAMING to " << newname << G4endl;
   renameViewPt(newname);
   listitem->setText(QString(newname));
   AuxWindowDialog->lineEdit->setText(newname);
   //   if (currentState == VIEWPOINT)
   //      scheduleRedraw();

   delete[] newname;
}

// Renames currently selected viewpoint.

void G4OpenInventorQtExaminerViewer::renameViewPt(char *vpName)
{
   std::size_t idx = 0, end, pos;
   std::size_t size = viewPtList.size();
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


void G4OpenInventorQtExaminerViewer::SortBookmarksCB()
{
   // FWJ NOTE: Xt version of this does not work (does nothing)

   //   G4cout << "Sort Button: SortBookmarks CALLBACK" << G4endl;

   // FWJ List for sorting
   // The dialog list and bookmark file will be rewritten.
   // Simpler to populate this list from the data structure.

   std::vector<std::string> charList;

   if (viewPtList.size() < 2) return;

   // Get current entries from the list

   for (int i = 0; i < (int)viewPtList.size(); i++) {

      charList.push_back(viewPtList[i].viewPtName);
      //      G4cout << "  Pushed " << i << " " << charList[i] << G4endl;
   }

   std::sort(charList.begin(), charList.end());

   // FWJ POPULATE the new dialog list
   //   G4cout << "  Populating Bookmark listWidget..." << G4endl;
   AuxWindowDialog->listWidget->clear();

   for (int i = 0; i < (int)viewPtList.size(); i++) {
      // viewPtIdx has to be changed to account for a different order in viewPtList
      if (!strcmp(charList[i].c_str(), curViewPtName))
         viewPtIdx = i;
      new QListWidgetItem(charList[i].c_str(), AuxWindowDialog->listWidget); 

   }

   sortViewPts(charList);

}

// Rewrites entire viewpoint file with sorted viewpoints.

void G4OpenInventorQtExaminerViewer::sortViewPts(std::vector<std::string> sortedViewPts) 
{
   SbVec3f axis;
   float x, y, z, angle;
   std::size_t sortIdx = 0, unsortIdx = 0;

   if (fileOut.is_open())
      fileOut.close();

   fileOut.open(fileName.c_str()); // Erase current viewpoint file

   writeViewPtIdx();

   std::size_t size = sortedViewPts.size();
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

// Needed to implement mouse wheel zoom direction change.
// Does not work with MacOS trackpad: use Coin3d default handler.
// Emulating private method SoGuiFullViewerP::zoom()
#ifndef __APPLE__
void
G4OpenInventorQtExaminerViewer::zoom(const float diffvalue)
{
   float multiplicator = float(std::exp(diffvalue));
   SoCamera *cam = getCamera();

   if (cam->isOfType(SoPerspectiveCamera::getClassTypeId())) {
      const float oldfocaldist = cam->focalDistance.getValue();
      const float newfocaldist = oldfocaldist * multiplicator;

      SbVec3f direction;
      cam->orientation.getValue().multVec(SbVec3f(0, 0, -1), direction);

      const SbVec3f oldpos = cam->position.getValue();
      const SbVec3f newpos = oldpos + (newfocaldist - oldfocaldist) * -direction;
      cam->position = newpos;
      cam->focalDistance = newfocaldist;
   } else if (cam->isOfType(SoOrthographicCamera::getClassTypeId())) {
      SoOrthographicCamera * oc = (SoOrthographicCamera *)cam;
      oc->height = oc->height.getValue() * multiplicator;
   }
}
#endif

// Handling mouse and keyboard events

SbBool
G4OpenInventorQtExaminerViewer::processSoEvent(const SoEvent* const ev)
{

   // FWJ DEBUG
   //   G4cout << "processSoEvent ############" << ++processSoEventCount << G4endl;

   SoCamera *cam = getCamera();
   const SoType type(ev->getTypeId());

// Needed to implement mouse wheel zoom direction change.
// Does not work with MacOS trackpad: use Coin3d default handler.
#ifndef __APPLE__
   if (type.isDerivedFrom(SoMouseButtonEvent::getClassTypeId())) {
      SoMouseButtonEvent * me = (SoMouseButtonEvent *) ev;

      //      if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
      //          || currentState == PAUSED_ANIMATION) {

      switch (me->getButton()) {

         case SoMouseButtonEvent::BUTTON4: // Scroll wheel up
            if (me->getState() == SoButtonEvent::DOWN) {
               //               G4cout << "SCROLL WHEEL UP" << G4endl;
               zoom(-0.1f);
               return TRUE;
            }
            break;

         case SoMouseButtonEvent::BUTTON5: // Scroll wheel down
            if (me->getState() == SoButtonEvent::DOWN) {
               //               G4cout << "SCROLL WHEEL DOWN" << G4endl;
               zoom(0.1f);
               return TRUE;
            }
            break;

         default:
            break;
      }
         //      }
      if (currentState == GENERAL) {

      }
   }
#endif

   if (type.isDerivedFrom(SoKeyboardEvent::getClassTypeId())) {
      SoKeyboardEvent* ke = (SoKeyboardEvent*)ev;

      if (SoKeyboardEvent::isKeyPressEvent(ev, ke->getKey())) {
         switch (ke->getKey()) {
         case SoKeyboardEvent::E:
            if (externalQtApp) {
               // G4cout << "E KEY PRESSED" << G4endl;
               return TRUE;
            } else {
               G4cout <<
                  "E KEY PRESSED, EXITING OIQT VIEWER SECONDARY LOOP" <<
                  G4endl;
               SoQt::exitMainLoop();
            //            escapeCallback();
               return TRUE;
            }
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
                    && currentState != REVERSED_ANIMATION
                    && currentState != PAUSED_ANIMATION)) {
               // FWJ IMPLEMENT LATER
               // String dialogName = (char *) "Delete Viewpoint";
               // std::string msg = "Are you sure you want to delete current viewpoint?";
               // warningMsgDialog(msg, dialogName, deleteViewPtCB);
               return TRUE;
            }
            break;
         case SoKeyboardEvent::LEFT_ARROW:
            switch (currentState) {
            case BEAMLINE:
               if ((this->lshiftdown)   || (this->rshiftdown)) {
                  refParticleIdx -= step;
                  moveCamera();
               }
               else if ((this->lctrldown)       || (this->rctrldown)) {
                  if (SoQtExaminerViewer::isAnimating())
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
               else {
                  if (SoQtExaminerViewer::isAnimating())
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
               SoDebugError::post("G4OpenInventorQtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::RIGHT_ARROW:
            switch(currentState) {
            case BEAMLINE:
               if ((this->lshiftdown)   || (this->rshiftdown)) {
                  refParticleIdx += step;
                  moveCamera();
               }
               else if ((this->lctrldown)       || (this->rctrldown)) {
                  if (SoQtExaminerViewer::isAnimating())
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
                  if (SoQtExaminerViewer::isAnimating())
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
               SoDebugError::post("G4OpenInventorQtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::DOWN_ARROW:
            switch(currentState) {
            case BEAMLINE:

               if ((this->lshiftdown)   || (this->rshiftdown)) {
                  refParticleIdx -= step;
                  moveCamera();
               }
               else{
                  if (SoQtExaminerViewer::isAnimating())
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
               SoDebugError::post("G4OpenInventorQtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::UP_ARROW:
            switch(currentState) {
            case BEAMLINE:
               if ((this->lshiftdown)   || (this->rshiftdown)) {
                  refParticleIdx -= step;
                  moveCamera();
               }
               else{
                  if (SoQtExaminerViewer::isAnimating())
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
               SoDebugError::post("G4OpenInventorQtExaminerViewer::processSoEvent",
                                  "Unhandled viewer state");
               break;
            }
            break;

         case SoKeyboardEvent::PAGE_UP:
            switch(currentState) {
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
               if(!animateSensor->isScheduled()) {
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
            default:    //fall through
               break;
            }
            break;

         case SoKeyboardEvent::PAGE_DOWN:
            switch(currentState) {
            case BEAMLINE:
               if (step > 1)
                  step--;
               return TRUE;
            case ANIMATION:
               if(!animateSensor->isScheduled()) {
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

            // FROM XT VIEWER
            //         case SoKeyboardEvent::E:
            //            this->escapeCallback(this->examinerObject);
            //            break;

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

   // Pass the event on to the viewer
   // Need some checks here as in Xt viewer?

   if (currentState == ANIMATION || currentState == REVERSED_ANIMATION
       || currentState == ROTATING)
      return FALSE;
   else
      return SoQtExaminerViewer::processSoEvent(ev);
}


// REMAINDER OF MENU BAR FUNCTIONS...


void G4OpenInventorQtExaminerViewer::FileLoadSceneGraphCB()
{
   //   G4cout << "File: Load scene graph CALLBACK" << G4endl;

   QFileDialog filedialog(getParentWidget(), tr("Load Scene Graph"));
   filedialog.setFileMode(QFileDialog::AnyFile);
   filedialog.setFont(*font);
   if (!filedialog.exec()) return;
   QStringList filenameinlist = filedialog.selectedFiles();
   QString filenamein = filenameinlist[0];

   //   G4cout << "Entered file name is " << qPrintable(filenamein) << G4endl;

   char* filename = new char[filenamein.size()+1];
   filename = strdup(qPrintable(filenamein));
   //   G4cout << "char[] file name is " << filename << G4endl;

   SoInput sceneInput;

   if (sceneInput.openFile(filename)) {
      // Read the whole file into the database
      newSceneGraph = SoDB::readAll(&sceneInput);
      if (newSceneGraph == NULL) {
         QMessageBox msgbox;
         msgbox.setFont(*font);
         QString messagetxt = "Error reading scene graph file ";
         messagetxt.append(filenamein);
         msgbox.setText(messagetxt);
         msgbox.exec();
         sceneInput.closeFile();
         return;
      }
   } else {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error opening scene graph file ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
      return;
   }

   SoSeparator* root = (SoSeparator*)getSceneGraph();
   root->unref();
   newSceneGraph->ref();
   setSceneGraph(newSceneGraph);
}

void G4OpenInventorQtExaminerViewer::FileSaveSceneGraphCB()
{
   //   G4cout << "File: Save scene graph CALLBACK" << G4endl;

   QFileDialog filedialog(getParentWidget(), tr("Save scene graph"));
   filedialog.setFileMode(QFileDialog::AnyFile);
   // To enable confirmation of overwriting
   filedialog.setAcceptMode(QFileDialog::AcceptSave);
   filedialog.setFont(*font);
   if (!filedialog.exec()) return;
   QStringList filenameinlist = filedialog.selectedFiles();
   QString filenamein = filenameinlist[0];

   //   G4cout << "Entered file name is " << qPrintable(filenamein) << G4endl;

   char* filename = new char[filenamein.size()+1];
   filename = strdup(qPrintable(filenamein));
   //   G4cout << "char[] file name is " << filename << G4endl;

   SoWriteAction writeAction;
   SoSeparator* root = (SoSeparator*)getSceneGraph();

   SoOutput* out = writeAction.getOutput();

   if (out->openFile(filename)) {
      out->setBinary(FALSE);
      writeAction.apply(root);
      out->closeFile();
   } else {
      QMessageBox msgbox;
      msgbox.setFont(*font);
      QString messagetxt = "Error opening file ";
      messagetxt.append(filenamein);
      msgbox.setText(messagetxt);
      msgbox.exec();
   }
}


void G4OpenInventorQtExaminerViewer::HelpControlsCB()
{
   //   G4cout << "Help: Help Controls CALLBACK" << G4endl;
   helpmsgbox->show();
}


HookEventProcState::HookEventProcState(G4OpenInventorQtExaminerViewer* vwr)
{
   viewer = vwr;
}

HookEventProcState::~HookEventProcState()
{;}

G4bool HookEventProcState::Notify(G4ApplicationState requestedState)
{
   if (requestedState == G4State_EventProc) viewer->newEvents = true;
   return true;
}
