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
//
#ifndef G4UIQt_h
#define G4UIQt_h 

#include <map>

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"

#include <qobject.h>
#include <qmap.h>
#include <qstringlist.h>
#include <qtabwidget.h>
#include <qdockwidget.h>
#include <qdialog.h>

class QMainWindow;
class QLineEdit;
class G4UIsession;
class QListWidget;
class QTreeWidget;
class QTreeWidgetItem;
class QTextEdit;
class QTextBrowser;
class QLabel;
class QResizeEvent;
class QTabWidget;
class QStringList;
class QSplitter;
class QToolBar;
class QTableWidget;
class QPixmap;
class QComboBox;
class QCompleter;
class QtGlobal;
class QStandardItemModel;
class QToolButton;

// Class description :
//
//  G4UIQt : class to handle a Qt interactive session.
// G4UIQt is the Qt version of G4UIterminal.
//
//  A command box is at disposal for entering/recalling Geant4 commands.
//  A menubar could be customized through the AddMenu, AddButton, AddIcon methods.
//  Note that there are corresponding Geant4 commands to add a 
// menus in the menubar and add buttons in a menu.
//  Ex : 
//    /gui/addMenu   test Test
//    /gui/addButton test Init /run/initialize
//    /gui/addButton test "Set gun" "/control/execute gun.g4m"
//    /gui/addButton test "Run one event" "/run/beamOn 1"
//
//  Command completion, by typing "tab" key, is available on the 
// command line.
//
// Class description - end :

class G4QTabWidget : public QTabWidget {
public :
  G4QTabWidget();
  G4QTabWidget(QWidget* aParent, G4int sizeX, G4int sizeY);
  void paintEvent  ( QPaintEvent * event );
  inline void setTabSelected(G4bool a) { fTabSelected = a; };
  inline void setLastTabCreated(G4int a) { fLastCreated = a; };
  inline bool isTabSelected() { return fTabSelected; };
  G4bool fTabSelected;
  G4int fLastCreated;
  G4int fPreferedSizeX;
  G4int fPreferedSizeY;
  inline void setPreferredSize(QSize s) {
    fPreferedSizeX = s.width() + 6; // tab label height + margin left+right
    fPreferedSizeY = s.height() + 58; // margin left+right
  }
  inline QSize sizeHint () const {
    return QSize(fPreferedSizeX, fPreferedSizeY);
  }
};

class G4UIOutputString {
  public :
  G4UIOutputString(QString text,G4String thread = "",G4String outputstream= "info");
  inline QString GetOutputList() { return " all info warning error ";};
  QString fText;
  G4String fThread;
  G4String fOutputStream; // Error, Warning, Info
};


class G4UIDockWidget : public QDockWidget {
public:
  G4UIDockWidget(QString txt);
  void closeEvent(QCloseEvent *);
};


class G4UIQt : public QObject, public G4VBasicShell, public G4VInteractiveSession {
  Q_OBJECT

public: // With description
  G4UIQt(G4int,char**);
  // (argv, argc) or (0, NULL) had to be given.
  G4UIsession* SessionStart();
  // To enter interactive X loop ; waiting/executing command,...
  void AddMenu(const char*,const char*);
  // To add a pulldown menu in the menu bar. 
  // First argument is the name of the menu.
  // Second argument is the label of the cascade button.
  // Ex : AddMenu("my_menu","My menu")
  void AddButton(const char*,const char*,const char*);
  // To add a push button in a pulldown menu.
  // First argument is the name of the menu.
  // Second argument is the label of the button.
  // Third argument is the Geant4 command executed when the button is fired.
  // Ex : AddButton("my_menu","Run","/run/beamOn 1"); 
  void AddIcon(const char* userLabel, const char* iconFile, const char* command, const char* file_name="");
  // To add a icon in the toolbar
  // First argument is the label of the icon.
  // Second argument is the selected icon type (open save move rotate pick zoom_in zoom_out wireframe solid hidden_line_removal hidden_line_and_surface_removal perspective ortho user_icon).
  // Third argument is the Geant4 command executed when the button is fired.
  // Fourth argument is the path to the icon file if "user_icon" selected
  // Ex : AddButton("change background color","../background.xpm"," /vis/viewer/set/background"); 
  void OutputStyle (const char*,const char*,const char*);
  // Specify an output style
  // First argument destination (cout cerr warnings errors all)
  // Second argument is the style (fixed proportional)
  // Third argument highlights commands if "highlight" (and if /control/verbose > 0)

  void NativeMenu(G4bool aVal);
  // Enable/Disable the native Menu Bar in Qt

  void ClearMenu();
  // Clear Menu Bar, remove all actions

  void DefaultIcons(G4bool aVal);
  // Enable/Disable the default icon ToolBar in Qt

  G4bool AddTabWidget(QWidget*,QString);
  // To add a tab for vis openGL Qt driver
  
  inline QTabWidget* GetViewerTabWidget() {
    return fViewerTabWidget;
  };

  QWidget* GetSceneTreeWidget();
  // Get the scene tree component

  QWidget* GetViewerPropertiesWidget();
  // Get the Viewer Properties Widget

  QWidget* GetPickInfosWidget();
  // Get the Pick Widget

  G4bool IsSplitterReleased();

  inline G4bool IsIconMoveSelected() {
    return fMoveSelected;
  };
  inline G4bool IsIconRotateSelected() {
    return fRotateSelected;
  };
  inline G4bool IsIconPickSelected() {
    return fPickSelected;
  };
  inline G4bool IsIconZoomInSelected() {
    return fZoomInSelected;
  };
  inline G4bool IsIconZoomOutSelected() {
    return fZoomOutSelected;
  };

  void SetIconMoveSelected();
  void SetIconRotateSelected();
  void SetIconPickSelected();
  void SetIconZoomInSelected();
  void SetIconZoomOutSelected();
  void SetIconHLHSRSelected();
  void SetIconHLRSelected();
  void SetIconSolidSelected();
  void SetIconWireframeSelected();
  void SetIconPerspectiveSelected();
  void SetIconOrthoSelected();

  inline QMainWindow * GetMainWindow() {
    return fMainWindow;
  };
  // Return the main window
  
  inline QPixmap* getSearchIcon() { return fSearchIcon;};
  // return the "search" icon pixmap
  inline QPixmap* getClearIcon() { return fClearIcon;};
  // return the "clear" icon pixmap

  void SetStartPage(const std::string&);
  // Set the text on the first page of the viewer. If "", will take the last value as default
  // Note: Qt Rich text format could be used, see link for example :
  // https://qt-project.org/doc/qt-4.8/richtext-html-subset.html#table-cell-attributes

  inline QWidget* GetCoutWidget() {
    return fCoutDockWidget->widget();
  };
  // Return the G4cout widget with filters

  inline G4UIDockWidget* GetCoutDockWidget() {
    return fCoutDockWidget;
  };
  // Return the cout dockable widget as a QDockWidget
  
  inline G4UIDockWidget* GetUserInterfaceWidget() {
    return fUIDockWidget;
  };
  // Return the UserInterface widget (including scene tree, help and History widgets)

  inline QTabWidget* GetUITabWidget() {
    return fUITabWidget;
  }
  // return the viewer widget including all viewers

  inline QWidget* GetHistoryWidget() {
    return fHistoryTBWidget;
  }
  // return the history widget
  
  inline QWidget* GetHelpWidget() {
    return fHelpTBWidget;
  }
  // return the help widget
  
  G4bool AddViewerTab(QWidget* w, std::string title);
  // Add a new tab in the viewer, could be used to add your own component
  
  G4bool AddViewerTabFromFile(std::string fileName, std::string title);
  // Add a new tab in the viewer containing the content of the file in a QLabel

public:
  ~G4UIQt();
  void Prompt(G4String);
  void SessionTerminate();
  virtual void PauseSessionStart(const G4String&);
  virtual G4int ReceiveG4cout(const G4String&);
  virtual G4int ReceiveG4cerr(const G4String&);
  //   G4String GetCommand(Widget);

private:
  void SecondaryLoop(G4String); // a VIRER
  void CreateHelpWidget();
  void InitHelpTreeAndVisParametersWidget();
  void FillHelpTree();
  void UpdateCommandCompleter();
  void CreateIcons();
  virtual void ExitHelp() const;
  void SetDefaultIconsToolbar();

  void CreateHelpTree(QTreeWidgetItem*,G4UIcommandTree*);
  QTreeWidgetItem* FindTreeItem(QTreeWidgetItem *,const QString&);

  QString GetCommandList(const G4UIcommand*);
  void updateHelpArea(const G4UIcommand*);
  virtual G4bool GetHelpChoice(G4int&);// have to be implemeted because we heritate from G4VBasicShell
  bool eventFilter(QObject*,QEvent*);
  void ActivateCommand(G4String);
#if QT_VERSION < 0x050F00
  QMap<G4int,QString> LookForHelpStringInChildTree(G4UIcommandTree *,const QString&);
#else
  QMultiMap<G4int,QString> LookForHelpStringInChildTree(G4UIcommandTree *,const QString&);
#endif
  QWidget* CreateVisParametersTBWidget();
  QWidget* CreateHelpTBWidget();
  G4UIDockWidget* CreateCoutTBWidget();
  QWidget* CreateHistoryTBWidget();
  G4UIDockWidget* CreateUITabWidget();
  QWidget* CreateSceneTreeWidget();
  void CreateViewerWidget();
  void OpenHelpTreeOnCommand(const QString &);
  QString GetShortCommandPath(QString);
  QString GetLongCommandPath(QTreeWidgetItem*);
  G4bool IsGUICommand(const G4UIcommand*);
  G4bool CreateVisCommandGroupAndToolBox(G4UIcommand*, QWidget*, G4int, G4bool isDialog);
  G4bool CreateCommandWidget(G4UIcommand* command, QWidget* parent, G4bool isDialog);
  void CreateViewerPropertiesDialog();
  void CreatePickInfosDialog();
#ifdef G4MULTITHREADED
  void UpdateCoutThreadFilter();
#endif
  void FilterAllOutputTextArea();
  QString FilterOutput(const G4UIOutputString&,const QString&,const QString&);
  G4String GetThreadPrefix();
  G4bool CheckG4EnvironmentVariable(char* txt, char* version);
  QStandardItemModel* CreateCompleterModel(G4String aCmd);
  void CreateEmptyViewerPropertiesWidget();
  void CreateEmptyPickInfosWidget();
private:

  QMainWindow * fMainWindow;
  QLabel *fCommandLabel;
  QLineEdit * fCommandArea;
  QTextEdit *fCoutTBTextArea;
  QTabWidget* fUITabWidget;
  std::vector <G4UIOutputString> fG4OutputString;
  QLineEdit * fCoutFilter;
  QCompleter* fCompleter;
  G4bool fDefaultIcons;
  
  QListWidget *fHistoryTBTableList;
  QTreeWidget *fHelpTreeWidget;
  QWidget* fHelpTBWidget;
  QWidget* fHistoryTBWidget;
  G4UIDockWidget* fCoutDockWidget;
  G4UIDockWidget* fUIDockWidget;
  QWidget* fSceneTreeWidget;
  QWidget* fViewerPropertiesWidget;
  QWidget* fPickInfosWidget;
  QLineEdit* fHelpLine;
  G4QTabWidget* fViewerTabWidget;
  QString fCoutText;
  QTextBrowser *fStartPage;
  QSplitter * fHelpVSplitter;
  QTextEdit* fParameterHelpLabel;
  QTableWidget* fParameterHelpTable;

  QToolBar *fToolbarApp;
  QToolBar *fToolbarUser;
  QString fStringSeparator;
  G4String fLastErrMessage;
  QString fLastOpenPath;
  QToolButton* fViewModePopupButton;
  QToolButton* fSurfaceModePopupButton;
  
  QPixmap* fSearchIcon;
  QPixmap* fClearIcon;
  QPixmap* fSaveIcon;
  QPixmap* fOpenIcon;
  QPixmap* fMoveIcon;
  QPixmap* fRotateIcon;
  QPixmap* fPickIcon;
  QPixmap* fZoomInIcon;
  QPixmap* fZoomOutIcon;
  QPixmap* fWireframeIcon;
  QPixmap* fSolidIcon;
  QPixmap* fHiddenLineRemovalIcon;
  QPixmap* fHiddenLineAndSurfaceRemovalIcon;
  QPixmap* fPerspectiveIcon;
  QPixmap* fOrthoIcon;
  QPixmap* fCommandIcon;
  QPixmap* fDirIcon;
  QPixmap* fRunIcon;
  QPixmap* fParamIcon;
  QPixmap* fPickTargetIcon;
  QPixmap* fExitIcon;
  
#ifdef G4MULTITHREADED
  QComboBox* fThreadsFilterComboBox;
#endif
  std::string fDefaultViewerFirstPageHTMLText;
  
  QDialog* fViewerPropertiesDialog;
  QDialog* fPickInfosDialog;
  QString fLastCompleteCommand;
  G4bool fMoveSelected;
  G4bool fRotateSelected;
  G4bool fPickSelected;
  G4bool fZoomInSelected;
  G4bool fZoomOutSelected;
  struct G4UIQtStyle {
    G4bool fixed, highlight;
  };
  std::map<G4String,G4UIQtStyle> fOutputStyles;

private Q_SLOTS :
  void ExitSession();
  void ClearButtonCallback();
  void SaveOutputCallback();
  void CommandEnteredCallback();
  void CommandEditedCallback(const QString & text);
  void ButtonCallback(const QString&);
  void HelpTreeClicCallback();
  void HelpTreeDoubleClicCallback();
  void ShowHelpCallback();
  void CommandHistoryCallback();
  void LookForHelpStringCallback();
  void UpdateTabWidget(G4int);
  void ResizeTabWidget( QResizeEvent* );
  void CoutFilterCallback(const QString&);
  void ThreadComboBoxCallback(G4int);
  void TabCloseCallback(G4int);
  void ToolBoxActivated(G4int);
  void VisParameterCallback(QWidget*);
  void ChangeColorCallback(QWidget*);
  void ChangeCursorAction(const QString&);
  void ChangeSurfaceStyle(const QString&);
  void OpenIconCallback(const QString&);
  void SaveIconCallback(const QString&);
  void ViewerPropertiesIconCallback(G4int);
  void ChangePerspectiveOrtho(const QString&);
};

#endif

