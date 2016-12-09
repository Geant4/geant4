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
// $Id: G4OpenGLQtViewer.hh 101105 2016-11-07 08:09:26Z gcosmo $
//
// 
// G4OpenGLQtViewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4
//
// 30/06/2014 : M.Kelsey : Change QPixmap objects to pointers

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLQTVIEWER_HH
#define G4OPENGLQTVIEWER_HH

#include "globals.hh"

#include "G4OpenGLViewer.hh"
#include "G4PhysicalVolumeModel.hh"

#include <qobject.h>
#include <qpoint.h>
#include <qpixmap.h>

class G4OpenGLSceneHandler;
class G4UImanager;
class G4Text;
class G4UIcommand;

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
class QSignalMapper;
class G4UIQt;
class QTableWidget;
class QTableWidgetItem;
class QScrollArea;
class QSplitter;

class G4OpenGLQtViewer: public QObject, virtual public G4OpenGLViewer {

  Q_OBJECT

    typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
    typedef std::vector<PVNodeID> PVPath;

public:
  G4OpenGLQtViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLQtViewer ();
#ifdef G4MULTITHREADED
  // In MT mode these functions are called in the following order for each run:
  // Called on the master thread before starting the vis sub-thread.
  virtual void DoneWithMasterThread ();
  // Called on the master thread after starting the vis sub-thread.
  virtual void MovingToVisSubThread ();
  // Called on the vis sub-thread when waiting for events.
  virtual void SwitchToVisSubThread ();
  // Called on the vis sub-thread when all events have been processed.
  virtual void DoneWithVisSubThread ();
  // Called on the vis sub-thread when all events have been processed.
  // virtual void MovingToMasterThread ();  Not used in G4OpenGLQtViewer.
  // Called on the master thread after the vis sub-thread has terminated.
  virtual void SwitchToMasterThread ();
#endif

private:
  G4OpenGLQtViewer (const G4OpenGLQtViewer&);
  G4OpenGLQtViewer& operator= (const G4OpenGLQtViewer&);
public:
  virtual void updateQWidget()=0;
  void updateSceneTreeWidget();
  void updateViewerPropertiesTableWidget();
  void updatePickInfosWidget(int, int);
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
  void addPVSceneTreeElement(const G4String& model,
                             G4PhysicalVolumeModel* pPVModel,
                             int currentPVPOIndex);
  void addNonPVSceneTreeElement(const G4String& model,
                                int currentPVPOIndex,
                                const std::string& modelDescription,
                                const G4Visible& visible);
  bool isTouchableVisible(int POindex);
  void clearTreeWidget();
  bool exportImage(std::string name="", int width=-1, int height=-1);

public:
  void G4MousePressEvent(QMouseEvent *event);
  void G4wheelEvent (QWheelEvent * event); 
  void G4keyPressEvent (QKeyEvent * event); 
  void G4keyReleaseEvent (QKeyEvent * event);
  void G4MouseDoubleClickEvent();
  void G4MouseReleaseEvent(QMouseEvent *evnt);
  void G4MouseMoveEvent(QMouseEvent *event);

protected:
  void CreateGLQtContext ();
  virtual void CreateMainWindow (QGLWidget*,const QString&);
  void G4manageContextMenuEvent(QContextMenuEvent *e);
  void rotateQtScene(float, float);
  void rotateQtSceneToggle(float, float);
  void moveScene(float, float, float,bool);
  void FinishView();
  void updateKeyModifierState(const Qt::KeyboardModifiers&);
  void displaySceneTreeComponent();
  G4Colour getColorForPoIndex(int poIndex);
  
  // So that privately accumulated vis attributes modifiers may be
  // concatenated with the standard vis attributes modifiers for commands
  // such as /vis/viewer/set/all and /vis/viewer/save...
  const std::vector<G4ModelingParameters::VisAttributesModifier>*
  GetPrivateVisAttributesModifiers() const;
  bool isCurrentWidget();

protected:
  QWidget* fGLWidget;
  bool hasPendingEvents();
  void savePPMToTemp();
  int fRecordFrameNumber;

  bool fHasToRepaint;
  bool fUpdateGLLock;
  bool fQGLWidgetInitialiseCompleted;
  bool fPaintEventLock;

  // Flag to indicate that action was initiated by interaction (mouse
  // click) on the scene tree.  It is used and reset in
  // G4OpenGLStoredQtViewer::CompareForKernelVisit to prevent rebuild
  // in this case.
  bool fMouseOnSceneTree;

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
  void setRecordingInfos(const QString&);
  QString getProcessErrorMsg();
  QWidget* getParentWidget();
  bool parseAndInsertInSceneTree(QTreeWidgetItem *,
                                  G4PhysicalVolumeModel* pPVModel,
                                 unsigned int fullPathIndex,
                                 const QString& parentRoot,
                                 unsigned int currentIndex,
                                 int currentPVPOIndex);
  void setCheckComponent(QTreeWidgetItem* item,bool check);
  void createSceneTreeComponent();
  void createSceneTreeWidget();
  void createViewerPropertiesWidget();
  void createPickInfosWidget();
  bool parseAndCheckVisibility(QTreeWidgetItem * treeNode,int POindex);
  QTreeWidgetItem* createTreeWidgetItem(const PVPath& fullPath,
                                     const QString& name,
                                     int copyNb,
                                     int POIndex,
                                     const QString& logicalName,
                                     Qt::CheckState state,
                                     QTreeWidgetItem * treeNode,
                                     const G4Colour& color);
  QString getModelShortName(const G4String& modelShortName);
  void cloneSceneTree(QTreeWidgetItem* rootItem);
  void changeDepthOnSceneTreeItem(double lookForDepth,double currentDepth,QTreeWidgetItem* item);
  void updatePositivePoIndexSceneTreeWidgetQuickMap(int POindex,QTreeWidgetItem* item);
  void changeQColorForTreeWidgetItem(QTreeWidgetItem* item, const QColor&);

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

// parse the scene tree and return a string of status that can be saved
  std::string parseSceneTreeAndSaveState();

  std::string parseSceneTreeElementAndSaveState(QTreeWidgetItem* item, unsigned int level);
  QString GetCommandParameterList (const G4UIcommand *aCommand);
  void changeColorAndTransparency(GLuint index, G4Color color);

#ifdef G4MULTITHREADED
  inline void SetQGLContextVisSubThread(QThread *th) {
    fQGLContextVisSubThread = th;
  }
  inline void SetQGLContextMainThread(QThread *th) {
    fQGLContextMainThread = th;
  }
#endif
  
  QMenu *fContextMenu;
  QPoint fLastPos1;
  QPoint fLastPos2;
  QPoint fLastPos3;
  QPoint fLastPickPoint;

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
  QAction *fMouseRotateAction;
  QAction *fMouseMoveAction;
  QAction *fMousePickAction;
  QAction *fMouseZoomInAction;
  QAction *fMouseZoomOutAction;
  QAction *fFullScreenOn;
  QAction *fFullScreenOff;
  QAction *fDrawingWireframe;
  QAction *fDrawingLineRemoval;
  QAction *fDrawingSurfaceRemoval;
  QAction *fDrawingLineSurfaceRemoval;
  QAction *fProjectionOrtho;
  QAction *fProjectionPerspective;
  G4OpenGLQtMovieDialog* fMovieParametersDialog;
  RECORDING_STEP fRecordingStep;
  QProcess *fProcess;
  QTime *fLastEventTime;
  int fSpinningDelay;
  int fNbMaxFramesPerSec;
  float fNbMaxAnglePerSec;
  int fLaunchSpinDelay;
  QWidget* fUISceneTreeWidget;
  QWidget* fUIViewerPropertiesWidget;
  QWidget* fUIPickInfosWidget;
  bool fNoKeyPress;
  bool fAltKeyPress;
  bool fControlKeyPress;
  bool fShiftKeyPress;
  bool fBatchMode;
  bool fCheckSceneTreeComponentSignalLock;
  bool fViewerPropertiesTableWidgetIsInit;
  QTreeWidget* fSceneTreeComponentTreeWidget;
  // This is only use to hold the old "expand" value, see file:///Developer/Documentation/Qt/html/qtreewidgetitem.html#setExpanded 
  QWidget* fSceneTreeWidget;
  bool fPVRootNodeCreate;
  QLineEdit* fFilterOutput;
  QString fFileSavePath;
  int fNbRotation ;
  int fTimeRotation;
  QString fTouchableVolumes;
  QDialog* fShortcutsDialog;
  QTableWidget *fViewerPropertiesTableWidget;
  QWidget* fPickInfosWidget;
  QScrollArea* fPickInfosScrollArea;
  int fTreeWidgetInfosIgnoredCommands;
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
  QSignalMapper *fSignalMapperMouse;
  QSignalMapper *fSignalMapperSurface;
  QSignalMapper *fSignalMapperPicking;

  // quick map index to find next item
  std::map <int, QTreeWidgetItem*>::const_iterator fLastSceneTreeWidgetAskForIterator;
  std::map <int, QTreeWidgetItem*>::const_iterator fLastSceneTreeWidgetAskForIteratorEnd;

  // quick map index to find next item
  std::map <int, QTreeWidgetItem*>::const_iterator fOldLastSceneTreeWidgetAskForIterator;
  std::map <int, QTreeWidgetItem*>::const_iterator fOldLastSceneTreeWidgetAskForIteratorEnd;

  // icons
  QPixmap* fTreeIconOpen;
  QPixmap* fTreeIconClosed;
  QPixmap* fSearchIcon;

  int fLastExportSliderValue;
  G4Color fLastHighlightColor;
  GLuint fLastHighlightName;
  bool fIsDeleting;

#ifdef G4MULTITHREADED
  QThread* fQGLContextVisSubThread;
  QThread* fQGLContextMainThread;
#endif
  
public Q_SLOTS :
  void startPauseVideo();

protected Q_SLOTS :
  void updateToolbarAndMouseContextMenu();

private Q_SLOTS :
  void actionSaveImage();
  void actionChangeBackgroundColor();
  void actionChangeTextColor();
  void actionChangeDefaultColor();
  void actionMovieParameters();

  void showShortcuts();
  void toggleMouseAction(int);
  void toggleSurfaceAction(int);
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
  void toggleSceneTreeComponentPickingCout(int);
  void togglePicking();
  void currentTabActivated(int);

  // action trigger by a click on a component scene tree
  void sceneTreeComponentSelected();
  void changeDepthInSceneTree(int);
  void changeSearchSelection();
  void changeColorAndTransparency(QTreeWidgetItem* item,int val);
  void tableWidgetViewerSetItemChanged(QTableWidgetItem *);
};

#endif

#endif
