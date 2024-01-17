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

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"
#include "G4SceneTreeItem.hh"

#include <qdialog.h>
#include <qdockwidget.h>
#include <qmap.h>
#include <qobject.h>
#include <qtabwidget.h>
#include <qtreewidget.h>

class QMainWindow;
class QLineEdit;
class G4UIsession;
class QListWidget;
class QTreeWidgetItem;
class QTextEdit;
class QTextBrowser;
class QLabel;
class QResizeEvent;
class QTabWidget;
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

class G4QTabWidget : public QTabWidget
{
 public:
  G4QTabWidget();
  G4QTabWidget(QWidget* aParent, G4int sizeX, G4int sizeY);
  void paintEvent(QPaintEvent* event) override;
  inline void setTabSelected(G4bool a) { fTabSelected = a; };
  inline void setLastTabCreated(G4int a) { fLastCreated = a; };
  inline bool isTabSelected() { return fTabSelected; };
  G4bool fTabSelected;
  G4int fLastCreated;
  G4int fPreferedSizeX;
  G4int fPreferedSizeY;
  inline void setPreferredSize(QSize s)
  {
    fPreferedSizeX = s.width() + 6;  // tab label height + margin left+right
    fPreferedSizeY = s.height() + 58;  // margin left+right
  }
  inline QSize sizeHint() const override { return QSize(fPreferedSizeX, fPreferedSizeY); }
};

class G4UIOutputString
{
 public:
  G4UIOutputString(QString text, G4String thread = "", G4String outputstream = "info");
  inline QString GetOutputList() { return " all info warning error "; };
  QString fText;
  G4String fThread;
  G4String fOutputStream;  // Error, Warning, Info
};

class G4UIDockWidget : public QDockWidget
{
 public:
  G4UIDockWidget(QString txt);
  void closeEvent(QCloseEvent*) override;
};

class G4UIQt : public QObject, public G4VBasicShell, public G4VInteractiveSession
{
  Q_OBJECT

 public:  // With description
  // (argv, argc) or (0, NULL) had to be given.
  G4UIQt(G4int, char**);

  // To enter interactive X loop ; waiting/executing command,...
  G4UIsession* SessionStart() override;

  // To add a pulldown menu in the menu bar.
  // First argument is the name of the menu.
  // Second argument is the label of the cascade button.
  // Ex : AddMenu("my_menu","My menu")
  void AddMenu(const char*, const char*) override;

  // To add a push button in a pulldown menu.
  // First argument is the name of the menu.
  // Second argument is the label of the button.
  // Third argument is the Geant4 command executed when the button is fired.
  // Ex : AddButton("my_menu","Run","/run/beamOn 1");
  void AddButton(const char*, const char*, const char*) override;

  // To add a icon in the toolbar
  // First argument is the label of the icon.
  // Second argument is the selected icon type (open save move rotate pick zoom_in zoom_out
  // wireframe solid hidden_line_removal hidden_line_and_surface_removal perspective ortho
  // user_icon). Third argument is the Geant4 command executed when the button is fired. Fourth
  // argument is the path to the icon file if "user_icon" selected Ex : AddButton("change background
  // color","../background.xpm"," /vis/viewer/set/background");
  void AddIcon(const char* userLabel, const char* iconFile, const char* command,
    const char* file_name = "") override;

  // Specify an output style - used by /gui/outputStyle
  // First argument destination ("cout" etc or "all")
  // Second argument is the required style - see guidance
  void SetOutputStyle(const char* destination, const char* style) override;

  // Enable/Disable the native Menu Bar in Qt
  void NativeMenu(G4bool aVal) override;

  // Clear Menu Bar, remove all actions
  void ClearMenu() override;

  // Enable/Disable the default icon ToolBar in Qt
  void DefaultIcons(G4bool aVal) override;

  // To add a tab for vis openGL Qt driver
  G4bool AddTabWidget(QWidget*, QString);

  inline QTabWidget* GetViewerTabWidget() { return fViewerTabWidget; };

  // Get the "old" scene tree component
  QWidget* GetSceneTreeWidget();

  // Get the Viewer Properties Widget
  QWidget* GetViewerPropertiesWidget();

  // Get the Pick Widget
  QWidget* GetPickInfosWidget();

  G4bool IsSplitterReleased();

  inline G4bool IsIconMoveSelected() { return fMoveSelected; };
  inline G4bool IsIconRotateSelected() { return fRotateSelected; };
  inline G4bool IsIconPickSelected() { return fPickSelected; };
  inline G4bool IsIconZoomInSelected() { return fZoomInSelected; };
  inline G4bool IsIconZoomOutSelected() { return fZoomOutSelected; };

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

  // Return the main window
  inline QMainWindow* GetMainWindow() { return fMainWindow; };

  // return the "search" icon pixmap
  inline QPixmap* getSearchIcon() { return fSearchIcon; };

  // return the "clear" icon pixmap
  inline QPixmap* getClearIcon() { return fClearIcon; };

  // Set the text on the first page of the viewer. If "", will take the last value as default
  // Note: Qt Rich text format could be used, see link for example :
  // https://qt-project.org/doc/qt-4.8/richtext-html-subset.html#table-cell-attributes
  void SetStartPage(const std::string&);

  // Return the G4cout widget with filters
  inline QWidget* GetCoutWidget() { return fCoutDockWidget->widget(); };

  // Return the cout dockable widget as a QDockWidget
  inline G4UIDockWidget* GetCoutDockWidget() { return fCoutDockWidget; };

  // Return the UserInterface widget (including scene tree, help and History widgets)
  inline G4UIDockWidget* GetUserInterfaceWidget() { return fUIDockWidget; };

  // return the viewer widget including all viewers
  inline QTabWidget* GetUITabWidget() { return fUITabWidget; }

  // return the history widget
  inline QWidget* GetHistoryWidget() { return fHistoryTBWidget; }

  // return the help widget
  inline QWidget* GetHelpWidget() { return fHelpTBWidget; }

  // Add a new tab in the viewer, could be used to add your own component
  G4bool AddViewerTab(QWidget* w, std::string title);

  // Add a new tab in the viewer containing the content of the file in a QLabel
  G4bool AddViewerTabFromFile(std::string fileName, std::string title);

  // Update "new" scene tree
  void UpdateSceneTree(const G4SceneTreeItem&) override;

public:
  ~G4UIQt() override;
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(const G4String&) override;
  G4int ReceiveG4debug(const G4String&) override;
  G4int ReceiveG4cout(const G4String&) override;
  G4int ReceiveG4cerr(const G4String&) override;
  //   G4String GetCommand(Widget);

 private:
  void SecondaryLoop(G4String);  // a VIRER
  void CreateHelpWidget();
  void InitHelpTreeAndVisParametersWidget();
  void FillHelpTree();
  void UpdateCommandCompleter();
  void CreateIcons();
  void ExitHelp() const override;
  void SetDefaultIconsToolbar();

  void CreateHelpTree(QTreeWidgetItem*, G4UIcommandTree*);
  QTreeWidgetItem* FindTreeItem(QTreeWidgetItem*, const QString&);

  // Create the "mother" widget
  QWidget* CreateSceneTreeWidget();

  // Classes/structs and functions for the "new" scene tree
  // UpdateSceneTree is in "public" section above.
  // Create and connect the new tree widget
  void CreateNewSceneTreeWidget();
  // Build Physical Volume tree of touchables
  void BuildPVQTree(const G4SceneTreeItem& g4stItem, QTreeWidgetItem* qtwItem);
  // Callbacks on new scene tree items
  void SceneTreeItemClicked(QTreeWidgetItem*);
  void SceneTreeItemDoubleClicked(QTreeWidgetItem*);
  void SceneTreeItemExpanded(QTreeWidgetItem*);
  void SceneTreeItemCollapsed(QTreeWidgetItem*);
  // Class for trapping special mouse events on new scene tree
  struct NewSceneTreeItemTreeWidget: public QTreeWidget {
    void mousePressEvent(QMouseEvent*) override;
    void ActWithoutParameter(const G4String& action, G4SceneTreeItem*);
    void ActWithABool(const G4String& action, G4SceneTreeItem*, G4bool);
    void ActWithAnInteger(const G4String& action, G4SceneTreeItem*);
    void ActWithADouble(const G4String& action, G4SceneTreeItem*);
    void ActWithAString(const G4String& action, G4SceneTreeItem*);
  };

  QString GetCommandList(const G4UIcommand*);
  void updateHelpArea(const G4UIcommand*);
  G4bool GetHelpChoice(
    G4int&) override;  // have to be implemeted because we heritate from G4VBasicShell
  bool eventFilter(QObject*, QEvent*) override;
  void ActivateCommand(G4String);
#if (QT_VERSION < QT_VERSION_CHECK(5, 15, 0))
  QMap<G4int, QString> LookForHelpStringInChildTree(G4UIcommandTree*, const QString&);
#else
  QMultiMap<G4int, QString> LookForHelpStringInChildTree(G4UIcommandTree*, const QString&);
#endif
  QWidget* CreateVisParametersTBWidget();
  QWidget* CreateHelpTBWidget();
  G4UIDockWidget* CreateCoutTBWidget();
  QWidget* CreateHistoryTBWidget();
  G4UIDockWidget* CreateUITabWidget();
  void CreateViewerWidget();
  void OpenHelpTreeOnCommand(const QString&);
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
  QString FilterOutput(const G4UIOutputString&, const QString&, const QString&);
  G4String GetThreadPrefix();
  G4bool CheckG4EnvironmentVariable(char* txt, char* version);
  QStandardItemModel* CreateCompleterModel(G4String aCmd);
  void CreateEmptyViewerPropertiesWidget();
  void CreateEmptyPickInfosWidget();

 private:
  QMainWindow* fMainWindow;
  QLabel* fCommandLabel;
  QLineEdit* fCommandArea;
  QTextEdit* fCoutTBTextArea;
  QTabWidget* fUITabWidget;
  std::vector<G4UIOutputString> fG4OutputString;
  QLineEdit* fCoutFilter;
  QCompleter* fCompleter;
  G4bool fDefaultIcons;

  QListWidget* fHistoryTBTableList;
  QTreeWidget* fHelpTreeWidget;
  QWidget* fHelpTBWidget;
  QWidget* fHistoryTBWidget;
  G4UIDockWidget* fCoutDockWidget;
  G4UIDockWidget* fUIDockWidget;
  QWidget* fSceneTreeWidget;
  QWidget* fNewSceneTreeWidget;
  NewSceneTreeItemTreeWidget* fNewSceneTreeItemTreeWidget;
  QWidget* fViewerPropertiesWidget;
  QWidget* fPickInfosWidget;
  QLineEdit* fHelpLine;
  G4QTabWidget* fViewerTabWidget;
  QString fCoutText;
  QTextBrowser* fStartPage;
  QSplitter* fHelpVSplitter;
  QTextEdit* fParameterHelpLabel;
  QTableWidget* fParameterHelpTable;

  QToolBar* fToolbarApp;
  QToolBar* fToolbarUser;
  QString fStringSeparator;
  G4String fLastErrMessage;
  QString fLastOpenPath;

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

 private Q_SLOTS:
  void ExitSession();
  void ClearButtonCallback();
  void SaveOutputCallback();
  void CommandEnteredCallback();
  void CommandEditedCallback(const QString& text);
  void ButtonCallback(const QString&);
  void HelpTreeClicCallback();
  void HelpTreeDoubleClicCallback();
  void ShowHelpCallback();
  void CommandHistoryCallback();
  void LookForHelpStringCallback();
  void UpdateTabWidget(int);
  void ResizeTabWidget(QResizeEvent*);
  void CoutFilterCallback(const QString&);
  void ThreadComboBoxCallback(int);
  void TabCloseCallback(int);
  void ToolBoxActivated(int);
  void VisParameterCallback(QWidget*);
  void ChangeColorCallback(QWidget*);
  void ChangeCursorAction(const QString&);
  void ChangeSurfaceStyle(const QString&);
  void OpenIconCallback(const QString&);
  void SaveIconCallback(const QString&);
  void ViewerPropertiesIconCallback(int);
  void ChangePerspectiveOrtho(const QString&);
};

#endif
