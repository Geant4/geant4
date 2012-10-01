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
// $Id: G4OpenGLQtViewer.hh,v 1.25 2010-10-08 10:07:31 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4OpenGLQtViewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLQTVIEWER_HH
#define G4OPENGLQTVIEWER_HH

#include "globals.hh"

#include "G4OpenGLViewer.hh"
#include "G4PhysicalVolumeModel.hh"

#include <qobject.h>
#include <qpoint.h>

class G4OpenGLSceneHandler;
class G4UImanager;
class G4Text;

class QGLWidget;
class QDialog;
class QTextEdit;
class QContextMenuEvent;
class QMenu;
class QImage;
class QAction;
class QTabWidget;
class QMouseEvent;
class QKeyEvent;
class QWheelEvent;
class QProcess;
class QTime;
class QVBoxLayout;
class QPushButton;
class QSlider;
class QTreeWidgetItem;
class QTreeWidget;
class QColor;
class G4OpenGLSceneHandler;
class G4OpenGLQtMovieDialog;
class QLineEdit;
class G4UIQt;

class G4OpenGLQtViewer: public QObject, virtual public G4OpenGLViewer {

  Q_OBJECT

    typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
    typedef std::vector<PVNodeID> PVPath;

public:
  G4OpenGLQtViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLQtViewer ();
private:
  G4OpenGLQtViewer (const G4OpenGLQtViewer&);
  G4OpenGLQtViewer& operator= (const G4OpenGLQtViewer&);
public:
  virtual void updateQWidget()=0;
  QString setEncoderPath(QString path);
  QString getEncoderPath();
  QString setTempFolderPath(QString path);
  QString getTempFolderPath();
  QString setSaveFileName(QString path);
  QString getSaveFileName();
  bool isRecording();
  bool isStopped();
  bool isPaused();
  bool isEncoding();
  bool isWaiting();
  bool isFailed();
  void setWaiting();
  bool isBadEncoder();
  bool isBadOutput();
  bool isBadTmp();
  bool isSuccess();
  void setBadTmp();
  void setBadOutput();
  void setBadEncoder();
  bool isReadyToEncode();
  void resetRecording();
  void encodeVideo();
  void stopVideo();
  void saveVideo();
  bool generateMpegEncoderParameters();
  void displayRecordingStatus();
  void DrawText(const G4Text&);
  void ResetView ();
  void addPVSceneTreeElement(const G4String model,
                             G4PhysicalVolumeModel* pPVModel,
                             int currentPVPOIndex);
  void addNonPVSceneTreeElement(const G4String model,
                                int currentPVPOIndex,
                                std::string modelDescription,
                                const G4Visible& visible);
  bool isTouchableVisible(int POindex);
  void clearTreeWidget();
public:
  void G4MousePressEvent(QMouseEvent *event);
  void G4wheelEvent (QWheelEvent * event); 
  void G4keyPressEvent (QKeyEvent * event); 
  void G4MouseDoubleClickEvent();
  void G4MouseReleaseEvent();
  void G4MouseMoveEvent(QMouseEvent *event);

protected:
  void CreateGLQtContext ();
  virtual void CreateMainWindow (QGLWidget*,QString);
  void G4manageContextMenuEvent(QContextMenuEvent *e);
  void rotateQtScene(float, float);
  void rotateQtSceneToggle(float, float);
  void moveScene(float, float, float,bool);
  void FinishView();
  void updateKeyModifierState(Qt::KeyboardModifiers);
  void displaySceneTreeComponent();
  G4Colour getColorForPoIndex(int poIndex);

protected:
  QGLWidget* fWindow;
  QWidget* fGLWindow;
  bool hasPendingEvents();
  void savePPMToTemp();
  int fRecordFrameNumber;

  bool fHasToRepaint;
  bool fReadyToPaint;
  bool fIsRepainting;

private:
  enum RECORDING_STEP {WAIT,START,PAUSE,CONTINUE,STOP,READY_TO_ENCODE,ENCODING,FAILED,SUCCESS,BAD_ENCODER,BAD_OUTPUT,BAD_TMP,SAVE}; 

  void createPopupMenu();
  void createRadioAction(QAction *,QAction *, const std::string&,unsigned int a=1);
  void rescaleImage(int, int);
  bool printPDF(const std::string,int,QImage);  
  void showMovieParametersDialog();
  void initMovieParameters();
  QString createTempFolder();
  QString removeTempFolder();
  void setRecordingStatus(RECORDING_STEP);
  void setRecordingInfos(QString);
  QString getProcessErrorMsg();
  QWidget* getParentWidget();
  bool parseAndInsertInSceneTree(QTreeWidgetItem *,
                                  G4PhysicalVolumeModel* pPVModel,
                                 unsigned int fullPathIndex,
                                 QString parentRoot,
                                 unsigned int currentIndex,
                                 int currentPVPOIndex);
  void setCheckComponent(QTreeWidgetItem* item,bool check);
  void initSceneTreeComponent();
  bool parseAndCheckVisibility(QTreeWidgetItem * treeNode,int POindex);
  QTreeWidgetItem* createTreeWidgetItem(PVPath fullPath,
                                     QString name,
                                     int copyNb,
                                     int POIndex,
                                     QString logicalName,
                                     Qt::CheckState state,
                                     QTreeWidgetItem * treeNode,
                                     G4Colour color);
  QString getModelShortName(G4String modelShortName);
  void cloneSceneTree(QTreeWidgetItem* rootItem);
  void changeDepthOnSceneTreeItem(double lookForDepth,double currentDepth,QTreeWidgetItem* item);
  void updatePositivePoIndexSceneTreeWidgetQuickMap(int POindex,QTreeWidgetItem* item);
  void changeQColorForTreeWidgetItem(QTreeWidgetItem* item, QColor);

  bool isSameSceneTreeElement(QTreeWidgetItem* parentOldItem,QTreeWidgetItem* parentNewItem);
  void changeOpenCloseVisibleHiddenSelectedColorSceneTreeElement(QTreeWidgetItem* subItem);
  bool isPVVolume(QTreeWidgetItem* item);
  QTreeWidgetItem* cloneWidgetItem(QTreeWidgetItem* item);
  void clearSceneTreeSelection(QTreeWidgetItem*);
  void clearTreeWidgetElements(QTreeWidgetItem* item);

  // Get the tree wigdet item for POindex if exists
  QTreeWidgetItem* getTreeWidgetItem(int POindex);

  // Get the old tree wigdet item for POindex if exists
  QTreeWidgetItem* getOldTreeWidgetItem(int POindex);


  QMenu *fContextMenu;
  QPoint fLastPos1;
  QPoint fLastPos2;
  QPoint fLastPos3;

  // delta of depth move. This delta is put in % of the scene view
  G4double fDeltaDepth;
  // delta of zoom move. This delta is put in % of the scene view
  G4double fDeltaZoom;
  // To ensure key event are keep one by one
  bool fHoldKeyEvent;
  // To ensure move event are keep one by one
  bool fHoldMoveEvent;
  // To ensure rotate event are keep one by one
  bool fHoldRotateEvent;
  bool fAutoMove;
  QString fEncoderPath;
  QString fTempFolderPath;
  QString fMovieTempFolderPath;
  QString fSaveFileName;
  QString fParameterFileName;
  QAction *fRotateAction;
  QAction *fMoveAction;
  QAction *fPickAction;
  QAction *fFullScreenOn;
  QAction *fFullScreenOff;
  QAction *fDrawingWireframe;
  QAction *fDrawingLineRemoval;
  QAction *fDrawingSurfaceRemoval;
  QAction *fDrawingLineSurfaceRemoval;
  G4OpenGLQtMovieDialog* fMovieParametersDialog;
  RECORDING_STEP fRecordingStep;
  QProcess *fProcess;
  QTime *fLastEventTime;
  int fSpinningDelay;
  int fNbMaxFramesPerSec;
  float fNbMaxAnglePerSec;
  int fLaunchSpinDelay;
  QTabWidget* fUISceneTreeComponentsTBWidget;
  bool fNoKeyPress;
  bool fAltKeyPress;
  bool fControlKeyPress;
  bool fShiftKeyPress;
  bool fBatchMode;
  bool fCheckSceneTreeComponentSignalLock;
  QTreeWidget* fSceneTreeComponentTreeWidget;
  // This is only use to hold the old "expand" value, see file:///Developer/Documentation/Qt/html/qtreewidgetitem.html#setExpanded 
  QTreeWidget* fOldSceneTreeComponentTreeWidget;
  QWidget* fSceneTreeWidget;
  bool fPVRootNodeCreate;
  QLineEdit* fHelpLine;


  int fNbRotation ;
  int fTimeRotation;
  QString fTouchableVolumes;
  QDialog* fShortcutsDialog;
  QTextEdit *fTreeInfoDialogInfos;
  QPushButton * fSceneTreeButtonApply;
  QTextEdit *fShortcutsDialogInfos;
  QSlider* fSceneTreeDepthSlider;
  std::map <int, PVPath > fTreeItemModels;
  std::map <int, PVPath > fOldTreeItemModels;

  // quick scene tree map
  std::map <int, QTreeWidgetItem*> fPositivePoIndexSceneTreeWidgetQuickMap;
  // old scene tree map
  std::map <int, QTreeWidgetItem*> fOldPositivePoIndexSceneTreeWidgetQuickMap;
  std::vector <QTreeWidgetItem*> fOldNullPoIndexSceneTreeWidgetQuickVector;
  // old vis attr color map
  std::map <int, QColor> fOldVisAttrColorMap;

  unsigned int fSceneTreeDepth;
  QTreeWidgetItem* fModelShortNameItem;
  int fNumber;
  int fMaxPOindexInserted;
  G4UIQt* fUiQt;

  // quick map index to find next item
  std::map <int, QTreeWidgetItem*>::const_iterator fLastSceneTreeWidgetAskFor;

  // quick map index to find next item
  std::map <int, QTreeWidgetItem*>::const_iterator fOldLastSceneTreeWidgetAskFor;


public Q_SLOTS :
  void startPauseVideo();

private Q_SLOTS :
  void actionMouseRotate();
  void actionMouseMove();
  void actionMousePick();
  void actionDrawingWireframe();
  void actionDrawingLineRemoval();
  void actionDrawingSurfaceRemoval();
  void actionDrawingLineSurfaceRemoval();
  void actionSaveImage();
  void actionChangeBackgroundColor();
  void actionChangeTextColor();
  void actionChangeDefaultColor();
  void actionMovieParameters();

  void showShortcuts();
  void toggleDrawingAction(int);
  void toggleRepresentation(bool);
  void toggleProjection(bool);
  void toggleTransparency(bool);
  void toggleAntialiasing(bool);
  void toggleHaloing(bool);
  void toggleAux(bool);
  void toggleHiddenMarkers(bool);
  void toggleFullScreen(bool);
  void processEncodeFinished();
  void processLookForFinished();
  void processEncodeStdout();
  void sceneTreeComponentItemChanged(QTreeWidgetItem* item, int id);
  void sceneTreeComponentSelected();
  void changeDepthInSceneTree(int);
  void changeSearchSelection();
  void changeColorAndTransparency(QTreeWidgetItem* item,int val);
  // Only use for Qt>4.0
  //  void dialogClosed();
};

#endif

#endif
