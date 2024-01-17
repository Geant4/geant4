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
// L. Garnier

#include "G4UIQt.hh"

#include "G4MTcoutDestination.hh"
#include "G4Qt.hh"
#include "G4StateManager.hh"
#include "G4Types.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#include "G4UIcommandTree.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4SceneTreeItem.hh"
#include "G4AttCheck.hh"

#include <qapplication.h>
#include <qdialog.h>
#include <qevent.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qpushbutton.h>
#include <qscrollbar.h>
#include <qsplitter.h>
#include <qtextbrowser.h>
#include <qtextedit.h>
#include <qwidget.h>

#include <cstring>
#include <qboxlayout.h>
#include <qbuttongroup.h>
#include <qcolordialog.h>
#include <qcombobox.h>
#include <qcompleter.h>
#include <qfiledialog.h>
#include <qgroupbox.h>
#include <qheaderview.h>
#include <qlistwidget.h>
#include <qmainwindow.h>
#include <qmenu.h>
#include <qpainter.h>
#include <qradiobutton.h>
#include <qscrollarea.h>
#include <qstandarditemmodel.h>
#include <qstringlist.h>
#include <qtabbar.h>
#include <qtablewidget.h>
#include <qtabwidget.h>
#include <qtextstream.h>
#include <qtoolbar.h>
#include <qtoolbox.h>

#include <QInputDialog>

#include <set>
#include <map>
#include <cstdlib>

#ifndef G4GMAKE
#  include "moc_G4UIQt.cpp"
#endif

// Pourquoi Static et non  variables de classe ?
static G4bool exitSession = true;
static G4bool exitPause = true;

/**   Build a Qt window with a menubar, output area and promt area<br>
<pre>
   +-----------------------+
   |exit menu|             |
   |                       |
   | +-------------------+ |
   | |                   | |
   | |  Output area      | |
   | |                   | |
   | +-------------------+ |
   |      | clear |        |
   | +-------------------+ |
   | |  promt history    | |
   | +-------------------+ |
   | +-------------------+ |
   | |> promt area       | |
   | +-------------------+ |
   +-----------------------+
</pre>
*/
G4UIQt::G4UIQt(G4int argc, char** argv)
  : fMainWindow(nullptr),
    fCommandLabel(nullptr),
    fCommandArea(nullptr),
    fCoutTBTextArea(nullptr),
    fUITabWidget(nullptr),
    fCoutFilter(nullptr),
    fCompleter(nullptr),
    fDefaultIcons(true),
    fHistoryTBTableList(nullptr),
    fHelpTreeWidget(nullptr),
    fHelpTBWidget(nullptr),
    fHistoryTBWidget(nullptr),
    fCoutDockWidget(nullptr),
    fUIDockWidget(nullptr),
    fSceneTreeWidget(nullptr),
    fNewSceneTreeWidget(nullptr),
    fNewSceneTreeItemTreeWidget(nullptr),
    fViewerPropertiesWidget(nullptr),
    fPickInfosWidget(nullptr),
    fHelpLine(nullptr),
    fViewerTabWidget(nullptr),
    fCoutText("Output"),
    fStartPage(nullptr),
    fHelpVSplitter(nullptr),
    fParameterHelpLabel(nullptr),
    fParameterHelpTable(nullptr),
    fToolbarApp(nullptr),
    fToolbarUser(nullptr),
    fStringSeparator("__$$$@%%###__"),
    fLastOpenPath(""),
    fSearchIcon(nullptr),
    fClearIcon(nullptr),
    fSaveIcon(nullptr),
    fOpenIcon(nullptr),
    fMoveIcon(nullptr),
    fRotateIcon(nullptr),
    fPickIcon(nullptr),
    fZoomInIcon(nullptr),
    fZoomOutIcon(nullptr),
    fWireframeIcon(nullptr),
    fSolidIcon(nullptr),
    fHiddenLineRemovalIcon(nullptr),
    fHiddenLineAndSurfaceRemovalIcon(nullptr),
    fPerspectiveIcon(nullptr),
    fOrthoIcon(nullptr),
    fCommandIcon(nullptr),
    fDirIcon(nullptr),
    fRunIcon(nullptr),
    fParamIcon(nullptr),
    fPickTargetIcon(nullptr),
    fExitIcon(nullptr)
#ifdef G4MULTITHREADED
    ,
    fThreadsFilterComboBox(nullptr)
#endif
    ,
    fDefaultViewerFirstPageHTMLText(""),
    fViewerPropertiesDialog(nullptr),
    fPickInfosDialog(nullptr),
    fLastCompleteCommand(""),
    fMoveSelected(false),
    fRotateSelected(true),
    fPickSelected(false),
    fZoomInSelected(false),
    fZoomOutSelected(false)
{
  G4Qt* interactorManager = G4Qt::getInstance(argc, argv, (char*)"Qt");
  if ((QApplication*)interactorManager->GetMainInteractor() == nullptr) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();

    if (verbose >= 2) {
      G4cout << "G4UIQt : Unable to init Qt. Aborted" << G4endl;
    }
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI != nullptr) UI->SetSession(this);
  if (UI != nullptr) UI->SetG4UIWindow(this);

  // Check if already define in external app QMainWindow
  G4bool found = false;
  Q_FOREACH (QWidget* widget, QApplication::allWidgets()) {
    if ((!found) && (widget->inherits("QMainWindow"))) {
      found = true;
    }
  }

  if (found) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();

    if (verbose >= 2) {
      G4cout << "G4UIQt : Found an external App with a QMainWindow already defined. Aborted"
             << G4endl;
    }
    return;
  }
  CreateIcons();

  fMainWindow = new QMainWindow();
  fMainWindow->setAttribute(Qt::WA_DeleteOnClose);

  fMainWindow->setCorner(Qt::TopLeftCorner, Qt::LeftDockWidgetArea);
  fMainWindow->setCorner(Qt::TopRightCorner, Qt::RightDockWidgetArea);
  fMainWindow->setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
  fMainWindow->setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);

  CreateViewerWidget();
  fMainWindow->addDockWidget(Qt::LeftDockWidgetArea, CreateUITabWidget());
  fMainWindow->addDockWidget(Qt::BottomDockWidgetArea, CreateCoutTBWidget());

  // Create the new scene tree stuff
  fNewSceneTreeWidget = new QWidget;
  fNewSceneTreeWidget->setStyleSheet ("padding: 0px ");
  fNewSceneTreeWidget->setLayout(new QVBoxLayout);
  fNewSceneTreeWidget->layout()->setContentsMargins(5,5,5,5);
  fNewSceneTreeWidget->setWindowTitle("some name"/*QString(GetName().data())*/);
  // Add it to the "old" fSceneTreeWidget
  fSceneTreeWidget->layout()->addWidget(fNewSceneTreeWidget);
  CreateNewSceneTreeWidget();

  // add defaults icons
  SetDefaultIconsToolbar();

  if (UI != nullptr) UI->SetCoutDestination(this);  // TO KEEP

#ifdef G4MULTITHREADED
  // explicitly request that cout/cerr messages from threads are ALSO propagated to the master.
  masterG4coutDestination = this;
#endif

  fMainWindow->setWindowTitle(QFileInfo(QCoreApplication::applicationFilePath()).fileName());
  fMainWindow->move(QPoint(50, 50));

  // force the size at be correct at the beggining
  // because the widget is not realized yet, the size of the main window is not up to date. But
  // we need it in order to add some viewer inside
  fMainWindow->resize(fUIDockWidget->width() + fCoutDockWidget->width() + 20,
    fUIDockWidget->height() + fCoutDockWidget->height() + 20);

  // set last focus on command line
  fCommandArea->setFocus(Qt::TabFocusReason);

  // Allow QTextCursor to be called by another thread :
  // http://qt-project.org/doc/qt-4.8/qmetatype.html#qRegisterMetaType
  qRegisterMetaType<QTextCursor>("QTextCursor");

  // add some tips
  AddTabWidget(fStartPage, "Useful tips");

  // Set not visible until session start
  fMainWindow->setVisible(false);
}

G4UIQt::~G4UIQt()
{
  G4UImanager* UI = G4UImanager::GetUIpointer();  // TO KEEP
  if (UI != nullptr) {  // TO KEEP
    UI->SetSession(nullptr);  // TO KEEP
    UI->SetG4UIWindow(nullptr);
    UI->SetCoutDestination(nullptr);  // TO KEEP
#ifdef G4MULTITHREADED
    masterG4coutDestination = nullptr;  // set to cout when UI is deleted
#endif
  }
}

void G4UIQt::DefaultIcons(G4bool aVal)
{
  fDefaultIcons = aVal;

  if (! fMainWindow->isVisible()) {
    return;
  }

  if (fToolbarApp != nullptr) {
    if (aVal) {
      fToolbarApp->setVisible(true);
    }
    else {
      // Set not visible until session start
      fToolbarApp->setVisible(false);
    }
  }
}

void G4UIQt::SetDefaultIconsToolbar()
{
  if (fDefaultIcons) {
    if (fToolbarApp == nullptr) {
      fToolbarApp = new QToolBar();
      fToolbarApp->setIconSize(QSize(20, 20));
      fMainWindow->addToolBar(Qt::TopToolBarArea, fToolbarApp);
    }

    // Open/Save Icons
    AddIcon("Open macro file", "open", "/control/execute");
    AddIcon("Save viewer state", "save", "/vis/viewer/save");

    // View parameters
    fToolbarApp->addAction(QIcon(*fParamIcon), "Viewer properties", this,
      [this]() { this->ViewerPropertiesIconCallback(0); });

    // Cursors style icons
    AddIcon("Move", "move", "");
    AddIcon("Pick", "pick", "");
    AddIcon("Zoom out", "zoom_out", "");
    AddIcon("Zoom in", "zoom_in", "");
    AddIcon("Rotate", "rotate", "");

    // Surface Style icons
    AddIcon("Hidden line removal", "hidden_line_removal", "");
    AddIcon("Hidden line and hidden surface removal", "hidden_line_and_surface_removal", "");
    AddIcon("Surfaces", "solid", "");
    AddIcon("Wireframe", "wireframe", "");

    // Perspective/Ortho icons
    AddIcon("Perspective", "perspective", "");
    AddIcon("Orthographic", "ortho", "");
    AddIcon("Run beam on", "runBeamOn", "/run/beamOn 1");
    AddIcon("Exit Application", "exit", "exit");
  }
}

// clang-format off
void G4UIQt::CreateIcons()
{
  const char * const save[]={
    "32 32 24 1",
    "       c None",
    "+      c #000200",
    "@      c #141E43",
    "#      c #000C56",
    "$      c #494A47",
    "%      c #636662",
    "&      c #312F2A",
    "*      c #191B19",
    "=      c #002992",
    "-      c #003DFF",
    ";      c #041DA5",
    ">      c #A8A9A3",
    ",      c #FDFFFC",
    "'      c #DDE0DD",
    ")      c #818783",
    "!      c #C9CBC8",
    "~      c #0116C3",
    "{      c #C5C8FA",
    "]      c #6596FC",
    "^      c #A0B4F9",
    "/      c #0B2AFD",
    "(      c #799BE3",
    "_      c #5F4826",
    ":      c #D5D8D5",
    "                                ",
    "                                ",
    "   +++++++++++++++++++++++++    ",
    "  +@##+$%%%%%%%%%%%%%%%&*$%&+   ",
    "  +=-;@>,,''',,,,,,,',,)&!,)+   ",
    "  +;-~@>,,,,,,,,,,,,,,,>$!,)+   ",
    "  +=-~@>,,,,,{]]]]]^,,,>*&$&+   ",
    "  +=-~@>,,,,,'{^{^^{,,,>*#=#+   ",
    "  +=-~@>,,,,,,,,,,,,,,,>@~/=+   ",
    "  +=-~@>,,,{{{''''{',,,>@~-=+   ",
    "  +=-~@>,,'^]]]]]]({,,,>@~-=+   ",
    "  +=-~@>,,,{{{{{{{{{,,,>@~-=+   ",
    "  +=-~@>,,,,,'{^{{^{,,,>@~-=+   ",
    "  +=-~@>,,,,,]]]]]]],,,>@~-=+   ",
    "  +=-~*>,,,,,,,,,,,,,,,>@~-=+   ",
    "  +=-~@>,,,,,,,,,,,,,,,>@~-=+   ",
    "  +=-/=$%%%%%%%%%%%%%%%$=/-=+   ",
    "  +=---;###############;---=+   ",
    "  +=---////////////////----=+   ",
    "  +=----------------///----=+   ",
    "  +=---=@##############@#--=+   ",
    "  +=---@+++++++++++*%))_+~-=+   ",
    "  +=---#+++++++++++&:,,>@~-=+   ",
    "  +=---#+++++++++++$',,>@~-=+   ",
    "  +=---#+++++++++++&!,,>@~-=+   ",
    "  +=/--#+++++++++++&',,>@~-=+   ",
    "   @;--#+++++++++++$',,>@~-=+   ",
    "    @;;@+++++++++++*)!>%@=;#+   ",
    "     @++++++++++++++*&**++@++   ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fSaveIcon = new QPixmap(save);

  const char * const search[]  = {
    /* columns rows colors chars-per-pixel */
    "19 19 8 1",
    "  c #5C5C5C",
    ". c #7D7D7D",
    "X c #9B9B9B",
    "o c #C3C3C3",
    "O c None",
    "+ c #000000",
    "@ c #000000",
    "# c None",
    /* pixels */
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOo.  .oOOOOOO",
    "OOOOOOX      XOOOOO",
    "OOOOOo  XOOX  oOOOO",
    "OOOOO. XOOOOX .OOOO",
    "OOOOO  OOOOOO  OOOO",
    "OOOOO  OOOOOO  OOOO",
    "OOOOO. XOOOOo .OOOO",
    "OOOOOo  oOOo  oOOOO",
    "OOOOOOX       XOOOO",
    "OOOOOOOo.  .   XOOO",
    "OOOOOOOOOOOOO.  XOO",
    "OOOOOOOOOOOOOO. XOO",
    "OOOOOOOOOOOOOOOoOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO"
  };
  fSearchIcon = new QPixmap(search);

  const char * const clear[]  = {
    /* columns rows colors chars-per-pixel */
    "20 20 8 1",
    "  c #020202",
    ". c #202020",
    "X c #2C2C2C",
    "o c #797979",
    "O c None",
    "+ c #797979",
    "@ c #797979",
    "# c #797979",
    /* pixels */
    "OOOOOOOOOOOOOOOOOOOO",
    "OOOOOOOo    oOOOOOOO",
    "OOOOOXX      XXOOOOO",
    "OOOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOOO",
    "OOOO XXXXXXXXXX OOOO",
    "OOO XOOOOOOOOOO  OOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOooOooOooO OOOO",
    "OOOOXOOOOOOOOOO OOOO",
    "OOOOOooooooooooOOOOO",
    "OOOOOO........OOOOOO"
  };

  fClearIcon = new QPixmap(clear);


  const char * const open[]={
    "32 32 33 1",
    "       c None",
    "+      c #09091E",
    "@      c #191B18",
    "#      c #5F615F",
    "$      c #777977",
    "%      c #AEB1AF",
    "&      c #929491",
    "*      c #515250",
    "=      c #858784",
    "-      c #333533",
    ";      c #000100",
    ">      c #272926",
    ",      c #424341",
    "'      c #696C6A",
    ")      c #5F4927",
    "!      c #583D18",
    "~      c #6E6A5B",
    "{      c #47351D",
    "]      c #E0A554",
    "^      c #FFD67B",
    "/      c #EFB465",
    "(      c #FDBF6C",
    "_      c #FFCD76",
    ":      c #806238",
    "<      c #362611",
    "[      c #0B0D0A",
    "}      c #68471B",
    "|      c #523E22",
    "1      c #B78A51",
    "2      c #A17B44",
    "3      c #D6A45E",
    "4      c #C29354",
    "5      c #A1A3A0",
    "                                ",
    "                                ",
    "                     +@@@#      ",
    "                    $%   +&   * ",
    "                   #=     $  -; ",
    "                           %>;+ ",
    "                           ,;;+ ",
    "  &#$''#'                 >;;;+ ",
    " =)!)!!!!~                *#$'' ",
    " {]^/((_({-  %%%%%%%%%%%        ",
    " {(^_^^^^:<{{{{{{{{{{{{{[&      ",
    " {/_/(((((/]]]]]]]]]]]/]!#      ",
    " {/^(((((_^^^^^^^^^^^^^^:#      ",
    " {/^(((_^^____________^^}$      ",
    " {/^(((((/////////////((!#      ",
    " {/^/^_:<|||||||||||||||@@****1 ",
    " {/^/^(<[)||||||||||||||))!!}<; ",
    " {/^_(:|234444444444444444432)1 ",
    " {/_^/<)34444444444444444443},  ",
    " {/^(2{:41111111111111111142|5  ",
    " {3^3<:31111111111111111143}-   ",
    " {/^2<:31111111111111111441|'   ",
    " {_/<:41111111111111111143},    ",
    " {(4<:31111111111111111144!#    ",
    " )4))44111111111111111144},     ",
    " )2<:31111111111111111144{#     ",
    " @|:14444444444444444444}*      ",
    " ;@434444444444444444434<#      ",
    " ;[))))))))))))))))))))!~       ",
    " ++++++++++++++++++++++;%       ",
    "                                ",
    "                                "}
  ;
  fOpenIcon = new QPixmap(open);


  const char * const move[]={
    "32 32 16 1",
    "       c None",
    ".      c #F1F1F1",
    "+      c #939393",
    "@      c #282828",
    "#      c #787878",
    "$      c #000000",
    "%      c #CCCCCC",
    "&      c #1A1A1A",
    "*      c #0D0D0D",
    "=      c #5D5D5D",
    "-      c #AEAEAE",
    ";      c #BBBBBB",
    ">      c #C9C9C9",
    ",      c #D6D6D6",
    "'      c #FFFFFF",
    ")      c #999999",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "               ..               ",
    "               ++               ",
    "              .@@.              ",
    "              #$$#              ",
    "             %&$$*%             ",
    "             =$$$$=             ",
    "            -**$$**-            ",
    "            %;%&*>;%            ",
    "          -%   @&   %-          ",
    "        ,=*;   @&   ;*=,        ",
    "      .#*$$>        >$$*#.      ",
    "    ')&$$$$*@@    @@*$$$$&)'    ",
    "    ')&$$$$*@@    @@*$$$$&+'    ",
    "      .#*$$>        >$$*#.      ",
    "        ,=*;   @&   ;*=,        ",
    "          -%   @&   %-          ",
    "            %;%&*>>%            ",
    "            -**$$**-            ",
    "             =$$$$=             ",
    "             %&$$*%             ",
    "              #$$#              ",
    "              .@@.              ",
    "               ++               ",
    "               ..               ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fMoveIcon = new QPixmap(move);

  const char * const rotate[]={
    "32 32 27 1",
    "       c None",
    ".      c #003333",
    "+      c #000066",
    "@      c #1A1A1A",
    "#      c #003399",
    "$      c #3333CC",
    "%      c #000033",
    "&      c #353535",
    "*      c #434343",
    "=      c #336699",
    "-      c #3399FF",
    ";      c #003366",
    ">      c #5D5D5D",
    ",      c #282828",
    "'      c #3399CC",
    ")      c #333333",
    "!      c #3366CC",
    "~      c #333399",
    "{      c #505050",
    "]      c #666666",
    "^      c #333366",
    "/      c #0033CC",
    "(      c #3366FF",
    "_      c #336666",
    ":      c #787878",
    "<      c #868686",
    "[      c #6B6B6B",
    "                   .++@         ",
    "                  #$$%&*        ",
    "                 =--; *>,       ",
    "                 '-=  )>&       ",
    "                !-',  ,>*       ",
    "             !!=--=    >*       ",
    "            =------!!~@&)@      ",
    "             --------!*{{{*&,   ",
    "             -------=){*{{{>>{) ",
    "            ,!-----=  ){&  ,&{{@",
    "          ,*>!----=   &>&     )@",
    "         ){>)~---=    *])      @",
    "        @*>,  --!     ,&@       ",
    "        @{*   '!      ,-!=~^,@  ",
    "        @&    ==      {/(----!^ ",
    "         _           ]:;(----'  ",
    "         ==_         >{+(----~  ",
    "          !-!!======!!(((---!   ",
    "           ='--------------!    ",
    "             =!!!!'!!=; !-!     ",
    "                   &<*  !~      ",
    "              @.  *[*   ;       ",
    "               ;+)>*            ",
    "                 @@             ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fRotateIcon = new QPixmap(rotate);

  const char * const pick[]={
    /* columns rows colors chars-per-pixel */
    "20 20 12 1 ",
    "  c #050804",
    ". c #222321",
    "X c #3B3C3A",
    "o c #4C4E4B",
    "O c #616360",
    "+ c #747673",
    "@ c #8A8C89",
    "# c #9FA19E",
    "$ c #BABCB9",
    "% c #CED0CD",
    "& c #E4E6E3",
    "* c None",
    /* pixels */
    "*********oo*********",
    "*********oo*********",
    "******$O.  .O%******",
    "****&o .O..O  O*****",
    "***&X @**oo**@ X****",
    "***o $***oo***$ O***",
    "**% @**********@ %**",
    "**O.***********& +**",
    "**.O*****@@*****o.**",
    "oo .oo**@  #*&XX. oo",
    "oo .oo**@  #*&oo. oO",
    "**.O*****##*****oX**",
    "**O ***********& +**",
    "**% @****&&****+ &**",
    "***O $***Xo***# +***",
    "****X @&*Xo*&+ o****",
    "*****O  o..o  +*****",
    "******%+.  X+&******",
    "*********oo*********",
    "*********oO*********"
  };
  fPickIcon = new QPixmap(pick);

  const char * const zoom_in[]={
    "32 32 11 1",
    "       c None",
    ".      c #C9CBC8",
    "+      c #A8A9A3",
    "@      c #818783",
    "#      c #D5D8D5",
    "$      c #9BCCCC",
    "%      c #5FC7F4",
    "&      c #FDFFFC",
    "*      c #636662",
    "=      c #9599CE",
    "-      c #DDE0DD",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "          .++@@++.              ",
    "         +++..#.+++             ",
    "       .@+...++++#+@.           ",
    "       @$.%%+&&&@%..@           ",
    "      ++.%%%+&&&*%%.++          ",
    "     .+#%%%%+&&&*%%.#+          ",
    "     ++..%%%+&&&*%%%.++         ",
    "     +#.+++++&&&*++++.+         ",
    "     @.+&&&&&&&&&&&&&+@         ",
    "     @#+&&&&&&&&&&&&&+@         ",
    "     @.+&&&&&&&&&&&&&+.         ",
    "     +++@***+&&&****@+.         ",
    "     ....++++&&&*++++..         ",
    "      ++.===+&&&*%=.++          ",
    "       @..==+&&&*=..@#&         ",
    "       .@+#.+&&&@-+@@*@         ",
    "         +++.++++++ *+@*        ",
    "          .+@@@++.  @**+*       ",
    "                    .*@*+*      ",
    "                     .*@*+*     ",
    "                      +*@@*     ",
    "                       .**+     ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fZoomInIcon = new QPixmap(zoom_in);

  const char * const zoom_out[]={
    "32 32 11 1",
    "       c None",
    ".      c #C9CBC8",
    "+      c #A8A9A3",
    "@      c #818783",
    "#      c #D5D8D5",
    "$      c #5FC7F4",
    "%      c #9BCCCC",
    "&      c #FDFFFC",
    "*      c #636662",
    "=      c #9599CE",
    "-      c #DDE0DD",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "          .++@@++.              ",
    "         +++..#.+++             ",
    "       .@+..$$$$.#+@.           ",
    "       @%.$$$$$$$$..@           ",
    "      ++.$$$$$$$$$$.++          ",
    "     .+#$$$$$$$$$$$.#+          ",
    "     ++..$$$$$$$$$$$.++         ",
    "     +#.+++++++++++++.+         ",
    "     @.+&&&&&&&&&&&&&+@         ",
    "     @#+&&&&&&&&&&&&&+@         ",
    "     @.+&&&&&&&&&&&&&+.         ",
    "     +++@***********@+.         ",
    "     ....++++++++++++..         ",
    "      ++.===$$$$$$=.++          ",
    "       @..===$$$$=..@#&         ",
    "       .@+#.$$$..-+@@*@         ",
    "         +++#--.+++ *+@*        ",
    "          .+@@@++.  @**+*       ",
    "                    .*@*+*      ",
    "                     .*@*+*     ",
    "                      +*@@*     ",
    "                       .**+     ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fZoomOutIcon = new QPixmap(zoom_out);

  const char * const wireframe[]={
    "32 32 24 1",
    "       c None",
    "+      c #E4E4E4",
    "@      c #D5D5D5",
    "#      c #E1E1E1",
    "$      c #E7E7E7",
    "%      c #D8D8D8",
    "&      c #A7A7A7",
    "*      c #000000",
    "=      c #989898",
    "-      c #8A8A8A",
    ";      c #B5B5B5",
    ">      c #1B1B1B",
    ",      c #676767",
    "'      c #959595",
    ")      c #4A4A4A",
    "!      c #878787",
    "~      c #D3D3D3",
    "{      c #C4C4C4",
    "]      c #A4A4A4",
    "^      c #5B5B5B",
    "/      c #B3B3B3",
    "(      c #787878",
    "_      c #C7C7C7",
    ":      c #585858",
    "                                ",
    "                  +@@#          ",
    "          $%@@@@@&****=+        ",
    "        +&********&@-***;       ",
    "   +@@@&**&@@@@@@$  @*-&>&+     ",
    "  +*****&+          %*@ ,**'#   ",
    "  @***)!~           @*{&*****+  ",
    "  @*!]***&+        +-*^**'~!*@  ",
    "  @*~ +@&**&@@@@@@&****&+  ~*@  ",
    "  @*@    +&********&-*=    @*@  ",
    "  @*@      $%@-*-@$ @*@    @*@  ",
    "  @*@         @*@   %*%    @*@  ",
    "  @*@         %*%   %*%    @*@  ",
    "  @*@         %*%   %*%    @*@  ",
    "  @*@         %*%   %*%    @*@  ",
    "  @*@         %*%   %*%    @*@  ",
    "  @*@         %*%   %*%    @*@  ",
    "  @*@         @*@   %*%    @*@  ",
    "  @*@         =*-+  @*@    @*@  ",
    "  @*@    $%@@&****&@-*-+   @*@  ",
    "  @*@ $@&*****&@@&******&~~!*@  ",
    "  @*{/***&@@%$    $@-*-&*****+  ",
    "  @*)*)(-~          @*@ ~)**]   ",
    "  +*******&@@@@+    %*_+]**]    ",
    "   +@@@@@&******&@%+_*^**]#     ",
    "          $%@@@&****:**&+       ",
    "                +%@&**&         ",
    "                    ++          ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fWireframeIcon = new QPixmap(wireframe);

  const char * const solid[]={
    "32 32 33 1",
    "       c None",
    "+      c #C2DEDE",
    "@      c #B5D7DF",
    "#      c #ACD6E6",
    "$      c #60C0EC",
    "%      c #4EB7EE",
    "&      c #53B9ED",
    "*      c #82CEEA",
    "=      c #CFDDDA",
    "-      c #94C9E8",
    ";      c #0960FF",
    ">      c #0943FF",
    ",      c #0949FF",
    "'      c #3CB3F0",
    ")      c #71C7EB",
    "!      c #73CBE5",
    "~      c #D3DDDB",
    "{      c #C4DDDE",
    "]      c #B7D5DF",
    "^      c #2DACF5",
    "/      c #59C1ED",
    "(      c #5FC0ED",
    "_      c #85CEE9",
    ":      c #096BFF",
    "<      c #2AACF6",
    "[      c #5CBEEC",
    "}      c #7ACAE4",
    "|      c #73CAEB",
    "1      c #71C8E5",
    "2      c #D1DDDA",
    "3      c #CBDDD9",
    "4      c #67C1EB",
    "5      c #80CDEA",
    "                                ",
    "                                ",
    "          +@@@@@@#$%&*=         ",
    "        +-;>>>>>>>>>,')!~       ",
    "   {]@@-;>>>>>>>>>>>>^/(_=      ",
    "  {:>>>>>>>>>>>>>>>>><//[)!=    ",
    "  ]>>>>>>>>>>>>>>>>>><////[)}   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>><//////|   ",
    "  @>>>>>>>>>>>>>>>>>></////[1   ",
    "  @>>>>>>>>>>>>>>>>>><////[*2   ",
    "  {:>>>>>>>>>>>>>>>>><//[)12    ",
    "   +@@@@@-;>>>>>>>>>><[)13      ",
    "          {]@@@-;>>>,'*3        ",
    "                +@@#452         ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fSolidIcon = new QPixmap(solid);

  const char * const hidden_line_removal[]={
    "32 32 15 1",
    "       c None",
    "+      c #D5D5D5",
    "@      c #C7C7C7",
    "#      c #9C9C9C",
    "$      c #000000",
    "%      c #8E8E8E",
    "&      c #808080",
    "*      c #A9A9A9",
    "=      c #D8D8D8",
    "-      c #CACACA",
    ";      c #181818",
    ">      c #9F9F9F",
    ",      c #ACACAC",
    "'      c #B9B9B9",
    ")      c #555555",
    "                                ",
    "                  +@@+          ",
    "          +@@@@@@#$$$$%+        ",
    "        +#$$$$$$$$#@&$$$*       ",
    "   =-@@#$$#@@@@@-=  @$&#;>=     ",
    "  =$$$$$#+          -$@ *$$%+   ",
    "  -$&@-=            -$-  #$$$=  ",
    "  -$@               -$-   +&$-  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    @$@  ",
    "  @$@               @$@    #$=  ",
    "  -$&@@@-=          -$-  =>;,   ",
    "  =$$$$$$$#@@@-=    -$'+#$$,    ",
    "   =-@@@@#$$$$$$#@-+'$)$$#+     ",
    "          =-@@@#$$$$)$$#+       ",
    "                +@@#$$#         ",
    "                    ++          ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fHiddenLineRemovalIcon = new QPixmap(hidden_line_removal);

  const char * const hidden_line_and_surface_removal[]={
    "32 32 40 1",
    "       c None",
    "+      c #FFFFFF",
    "@      c #89A2E9",
    "#      c #5378E3",
    "$      c #A2B5ED",
    "%      c #5379E3",
    "&      c #5076E3",
    "*      c #3E69E4",
    "=      c #0C43F8",
    "-      c #043FFE",
    ";      c #CDD9ED",
    ">      c #BDCDE9",
    ",      c #FBFCFC",
    "'      c #406AE4",
    ")      c #0439FE",
    "!      c #0137FF",
    "~      c #4F75E3",
    "{      c #9EB5E3",
    "]      c #829FE0",
    "^      c #B6C6E7",
    "/      c #9DB4E3",
    "(      c #7E9CE0",
    "_      c #B2C3E9",
    ":      c #7E9AE0",
    "<      c #86A2E1",
    "[      c #CAD6ED",
    "}      c #5177E3",
    "|      c #829CE0",
    "1      c #BCCCE9",
    "2      c #3A67E6",
    "3      c #0A43FA",
    "4      c #95ACE1",
    "5      c #BBCBE9",
    "6      c #A9BBE5",
    "7      c #96AFE1",
    "8      c #BDCBE9",
    "9      c #4067E4",
    "0      c #6485E5",
    "a      c #E3EAF3",
    "b      c #CAD6F3",
    "                                ",
    "                                ",
    "                  ++++          ",
    "          ++++++++@#$+++        ",
    "        ++@%####&*=-#+;>,       ",
    "   +++++@'=)))))))!)~+{]^++     ",
    "   +$%&*=)!!!!!!!!!)~+/(]_+++   ",
    "   +#-))!!!!!!!!!!!)~+/(::<[+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/::::{+   ",
    "   +#)!!!!!!!!!!!!!!}+/:::|1+   ",
    "   +$#}}~23!!!!!!!!)~+/(]45,    ",
    "   +++++++@#}}~23!!)~+678++     ",
    "          ++++++@#~90+a++       ",
    "                ++++b++         ",
    "                    ++          ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fHiddenLineAndSurfaceRemovalIcon = new QPixmap(hidden_line_and_surface_removal);

  const char * const perspective[]={
    "32 32 3 1",
    "       c None",
    ".      c #D5D8D5",
    "+      c #000000",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "                                ",
    "           ................     ",
    "       ....+++++++++++++++.     ",
    "    ...++++..+.........+++.     ",
    "   ..++..............++..+.     ",
    "   .+++++++++++++++++.. .+.     ",
    "   .+...............+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+.    .+.  .+.     ",
    "   .+.      .+......+....+.     ",
    "   .+.     ..++++++.+.++++.     ",
    "   .+.    .++.......+...+..     ",
    "   .+.   .++.      .+..++.      ",
    "   .+. ..+..       .+..+.       ",
    "   .+..++.         .+.+.        ",
    "   .+.++.          .+++.        ",
    "   .+++.............++.         ",
    "   .+++++++++++++++++.          ",
    "   ...................          ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fPerspectiveIcon = new QPixmap(perspective);

  const char * const ortho[]={
    "32 32 3 1",
    "       c None",
    ".      c #D5D8D5",
    "@      c #000000",
    "                                ",
    "                                ",
    "                                ",
    "          ...................   ",
    "         ..@@@@@@@@@@@@@@@@@.   ",
    "       ..@@@.............@@@.   ",
    "      ..@@.@.         ..@..@.   ",
    "    ..@@ ..@.        .@@...@.   ",
    "   ..@@..............@@.. .@.   ",
    "   .@@@@@@@@@@@@@@@@@..   .@.   ",
    "   .@...............@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@.      .@.    .@.   ",
    "   .@.    .@........@......@.   ",
    "   .@.   .@@@@@@@@@.@.@@@@@@.   ",
    "   .@.  .@@+........@....@@..   ",
    "   .@...@.         .@...@...    ",
    "   .@.@@.          .@.@@ .      ",
    "   .@@@.............@@@..       ",
    "   .@@@@@@@@@@@@@@@@@...        ",
    "   ...................          ",
    "                                ",
    "                                ",
    "                                "}
  ;
  fOrthoIcon = new QPixmap(ortho);

  const char * const commandIcon[]={
    "20 20 25 1 ",
    "  c #4ED17F",
    ". c #4FD280",
    "X c #50D381",
    "o c #5BD181",
    "O c #5DD382",
    "+ c #59D48A",
    "@ c #66D68C",
    "# c #6FD895",
    "$ c #85DEA4",
    "% c #8CE0AC",
    "& c #96E4B8",
    "* c #9EE3B8",
    "= c #A8E5BB",
    "- c #A7E8C4",
    "; c #B2EAC8",
    ": c #B9ECD1",
    "> c #C2EDD3",
    ", c #CBF1DF",
    "< c #D4F3E3",
    "1 c #DDF4E5",
    "2 c #DBF5EC",
    "3 c #E5F7F0",
    "4 c #EDFAFB",
    "5 c #F6FBFE",
    "6 c #FEFFFC",
    /* pixels */
    "66666666666666666666",
    "66%++++++++++++++&56",
    "6$ o..o......o..o *6",
    "6+o...o*<441;@.o..+6",
    "6+..o@1553<354$..o+6",
    "6+..o<5<@  .*54#o.+6",
    "6+o.*52X     :5-..@6",
    "6+..15%      o$+o.+6",
    "6+.+55@        .o.+6",
    "6O.#54         .X.+6",
    "6O #54         .X.+6",
    "6O.+55@        .o.+6",
    "6+..25%      @,*o.@6",
    "6+o.*52X     :5>.o+6",
    "6+..O25<@  X=54#o.+6",
    "6+.o.@1553<354$...@6",
    "6+o..oo*<44<;@o..o+6",
    "6$ .o..o.....o..o *6",
    "66%+++++OOOO+++++*66",
    "66666666666666666666"
  };
  fCommandIcon = new QPixmap(commandIcon);

  const char * const dirIcon[]={
    "20 20 25 1 ",
    "  c #DF5959",
    ". c #DD5F5F",
    "X c #DE7370",
    "o c #E06360",
    "O c #E06467",
    "+ c #E06C6C",
    "@ c #E57979",
    "# c #E08886",
    "$ c #E18D91",
    "% c #E19D9B",
    "& c #E99B9D",
    "* c #E8A2A2",
    "= c #EEB2B0",
    "- c #EDBBBC",
    "; c #EDCBC7",
    ": c #E9CDD1",
    "> c #F1D5D6",
    ", c #F9DFE2",
    "< c #EFE8E7",
    "1 c #F3E3E4",
    "2 c #F8EEEC",
    "3 c #FCF6F4",
    "4 c #F6F3F9",
    "5 c #F2F8FC",
    "6 c #FEFFFD",
    /* pixels */
    "66666666666666666666",
    "66$oOOOOOOOOOOOOo%66",
    "6#                %6",
    "6o  +,666663:+    o6",
    "6o  =635533666$   o6",
    "6o  -65:+  +165X  o6",
    "6o  >6<.     36;  O6",
    "6o  26-      &6>. o6",
    "6. o56*      @63. o6",
    "6. X56&      o66. o6",
    "6. X56&      +63. o6",
    "6. o56*      @62. o6",
    "6o  26-      =61  O6",
    "6o  >6<.    o36:  o6",
    "6o  -65:+  @265X  o6",
    "6o  =635543665#   O6",
    "6o  +1666662;+    o6",
    "6#                %6",
    "66$OOOoo....OOOOo%66",
    "66666666666666666666"}
  ;
  fDirIcon = new QPixmap(dirIcon);


  const char * const runIcon[]={
    /* columns rows colors chars-per-pixel */
    "20 20 33 1 ",
    "  c #5CA323",
    ". c #5EA03F",
    "X c #6DB620",
    "o c #66AD3F",
    "O c #70B73C",
    "+ c #7CC13F",
    "@ c #569B41",
    "# c #61A14E",
    "$ c #70A95D",
    "% c #7EB55C",
    "& c #85B94E",
    "* c #90BE49",
    "= c #81B669",
    "- c #81B370",
    "; c #95CA46",
    ": c #A1CD40",
    "> c #AED045",
    ", c #B3D558",
    "< c #9BC87E",
    "1 c #AED668",
    "2 c #A2D075",
    "3 c #C2DC73",
    "4 c #A5C98F",
    "5 c #C1DC9F",
    "6 c #CAE18E",
    "7 c #CCE39A",
    "8 c #C4DCB6",
    "9 c #E3ECBA",
    "0 c #EEF3D3",
    "q c #F0F7DE",
    "w c #F8FAE9",
    "e c #FCFFFB",
    "r c None",
    /* pixels */
    "rrrrrrrr%<<2rrrrrrrr",
    "rrrrr5=$$$$===rrrrrr",
    "rrrr<##$$$$$---&rrrr",
    "rrr=###$$$$-----%rrr",
    "rr=####$$$$------&rr",
    "r2@####7##$-------rr",
    "r.@####048$-------Or",
    "r.@####q4ee=----$@.r",
    " .@@###w4eee5%$#@@@X",
    " .@@@..w4eeeeqo..@@X",
    " .@..ooe<eeee7Oooo@X",
    " ..oooOe2eee6OOOooo ",
    "rOooOO+e2ew2+++++O+r",
    "r:oO+++e30,;;;;;++Or",
    "r :++;:9,>,,>>:;;1rr",
    "rr*1;:>,333333,>32rr",
    "rrr66,1367777637<rrr",
    "rrrr509799999905rrrr",
    "rrrrr=8wqwwww8-rrrrr",
    "rrrrrrrr4444rrrrrrrr"
  };
  fRunIcon = new QPixmap(runIcon);

  const char * const paramIcon[]={
    /* columns rows colors chars-per-pixel */
    "20 20 35 1 ",
    "  c #2E2525",
    ". c #403737",
    "X c #423A3A",
    "o c #443C3C",
    "O c #473F3F",
    "+ c #4C4444",
    "@ c #4F4848",
    "# c #514949",
    "$ c #544D4D",
    "% c #595252",
    "& c #625B5B",
    "* c #696262",
    "= c #6D6666",
    "- c #716B6B",
    "; c #726C6C",
    ": c #767171",
    "> c #7E7878",
    ", c #8B8787",
    "< c #8C8787",
    "1 c #8D8888",
    "2 c #918D8D",
    "3 c #928E8E",
    "4 c #948F8F",
    "5 c #9C9898",
    "6 c #9D9999",
    "7 c #D5D4D4",
    "8 c #D8D6D6",
    "9 c #DDDBDB",
    "0 c #EFEFEF",
    "q c #F6F6F6",
    "w c None",
    "e c None",
    "r c None",
    "t c gray99",
    "y c None",
    /* pixels */
    "wwwwwwww5  5wwwwwwww",
    "wwwwwwww,  ,wwwwwwww",
    "www&;ww7+  +9ww=-www",
    "ww&  O#      OX  *ww",
    "ww;              >ww",
    "wwwO    .%%X    +www",
    "www#   3wwww3   Owww",
    "ww7   3wwwwww3   7ww",
    "5<+  .wwwwwww0.  +<5",
    "     %wwwwwwww$     ",
    "     %wwwwwwww$     ",
    "5<+  .wwwwwww0X  +<5",
    "ww9   4wwwwww1   9ww",
    "wwwO   30ww03   Owww",
    "wwwX    X#$X    @www",
    "ww=              =ww",
    "ww-  +O      ++  :ww",
    "www*>ww7+  +7ww=:www",
    "wwwwwwww1  1wwwwwwww",
    "wwwwwwww5  5wwwwwwww"
  };
  fParamIcon = new QPixmap(paramIcon);

  const char * const exitIcon[]={
    /* columns rows colors chars-per-pixel */
    "23 28 55 1 ",
    "  c None",
    ". c #350505",
    "X c #3A0505",
    "o c #3C0605",
    "O c #3D0605",
    "+ c #430606",
    "@ c #440606",
    "# c #470706",
    "$ c #500707",
    "% c #510807",
    "& c #520807",
    "* c #530807",
    "= c #550808",
    "- c #570808",
    "; c #5C0908",
    ": c #5D0908",
    "> c #5F0908",
    ", c #630A08",
    "< c #640A09",
    "1 c #6B0A09",
    "2 c #6C0A09",
    "3 c #720B0A",
    "4 c #760B0A",
    "5 c #770B0A",
    "6 c #7A0B0B",
    "7 c #7D0C0B",
    "8 c #7F0C0B",
    "9 c #840D0B",
    "0 c #850D0C",
    "q c #880D0C",
    "w c #8D0E0C",
    "e c #900E0C",
    "r c #940E0D",
    "t c #950E0D",
    "y c #9C0F0E",
    "u c #9E100E",
    "i c #AA100E",
    "p c #AC100F",
    "a c #AD100F",
    "s c #AE110F",
    "d c #B31110",
    "f c #B51210",
    "g c #B61210",
    "h c #B71210",
    "j c #B91210",
    "k c #C01311",
    "l c #C21311",
    "z c #C81311",
    "x c #C91312",
    "c c #CC1412",
    "v c #CE1412",
    "b c #D01412",
    "n c #D11412",
    "m c #D31412",
    "M c #D51513",
    /* pixels */
    "                       ",
    "          O=           ",
    "         :MMh          ",
    "         hMMM          ",
    "         jMMM          ",
    "    <x1  jMMM  %xw     ",
    "   rMMM  jMMM  MMMk    ",
    "  rMMMM# jMMM  MMMMx   ",
    " OMMMMk  jMMM  8MMMM9  ",
    " xMMMM   jMMM   pMMMM  ",
    " MMMM    jMMM    xMMM8 ",
    "rMMM3    jMMM     MMMM ",
    "MMMM     hMMM     MMMM ",
    "MMMM     :MMh     hMMM ",
    "MMMM      O%      8MMM ",
    "MMMM              pMMM ",
    "MMMM              MMMM ",
    "wMMM3             MMMM ",
    ".MMMM            xMMM9 ",
    " hMMMk          wMMMM  ",
    "  MMMMMO       hMMMM=  ",
    "  <MMMMMp:  $rMMMMMp   ",
    "   yMMMMMMMMMMMMMMk    ",
    "    #MMMMMMMMMMMM3     ",
    "      uMMMMMMMMk       ",
    "        #1wr3%         ",
    "                       ",
    "                       "
  };
  fExitIcon= new QPixmap(exitIcon);
}
// clang-format on

namespace {
  G4SceneTreeItem* ConvertToG4SceneTreeItem(QTreeWidgetItem* item)
  {
    auto qVariant = item->data(0, Qt::UserRole);
    std::istringstream iss(qVariant.toString().toStdString());
    void* itemAddress; iss >> itemAddress;
    return static_cast<G4SceneTreeItem*>(itemAddress);
  }

  QColor ConvertG4ColourToQColor(const G4Colour& g4Colour)
  {
    return QColor((int)(g4Colour.GetRed()*255),
                  (int)(g4Colour.GetGreen()*255),
                  (int)(g4Colour.GetBlue()*255),
                  (int)(g4Colour.GetAlpha()*255));
  }

  G4Colour ConvertQColorToG4Colour(const QColor& qColor)
  {
    return G4Color(qColor.red()/255.,
                   qColor.green()/255.,
                   qColor.blue()/255.,
                   qColor.alpha()/255.);
  }
}

void G4UIQt::CreateNewSceneTreeWidget()
{
  auto vLayout = fNewSceneTreeWidget->layout();
  // reduce margins
  vLayout->setContentsMargins(0,0,0,0);

  fNewSceneTreeItemTreeWidget = new NewSceneTreeItemTreeWidget;
  fNewSceneTreeItemTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  vLayout->addWidget(fNewSceneTreeItemTreeWidget);

  // A click on the item is handled here.
  // A click on the check box makes the volume visible/invisible
  connect(fNewSceneTreeItemTreeWidget, &QTreeWidget::itemClicked,
          [&](QTreeWidgetItem* item){SceneTreeItemClicked(item);});
  // A double click on either the colour icon or the name pops up a colour dialogue
  connect(fNewSceneTreeItemTreeWidget, &QTreeWidget::itemDoubleClicked,
          [&](QTreeWidgetItem* item){SceneTreeItemDoubleClicked(item);});

  // A click on the expansion triangle is handled here.
  connect(fNewSceneTreeItemTreeWidget, &QTreeWidget::itemExpanded,
          [&](QTreeWidgetItem* item){SceneTreeItemExpanded(item);});
  connect(fNewSceneTreeItemTreeWidget, &QTreeWidget::itemCollapsed,
          [&](QTreeWidgetItem* item){SceneTreeItemCollapsed(item);});
}
  
void G4UIQt::UpdateSceneTree(const G4SceneTreeItem& root)
{
  //  G4debug << "\nG4UIQt::UpdateSceneTree: scene tree summary\n";
  //  root.DumpTree(G4debug);  // Verbosity = 0 (one line per item)
  //  G4debug << "\nG4UIQt::UpdateSceneTree: scene tree dump\n";
  //  root.DumpTree(G4debug,1);  // Verbosity = 1 (higher levels available)

  // Clear the existing GUI-side tree
  fNewSceneTreeItemTreeWidget->clear();
  // (I think this deletes everything - the top level items and their children.)

  // Build a new GUI-side tree
  fNewSceneTreeItemTreeWidget->setHeaderLabel (root.GetDescription().c_str());
  for (const auto& model : root.GetChildren()) {

    auto item = new QTreeWidgetItem(fNewSceneTreeItemTreeWidget);

    // Add this GUI-side representation of the model as a child of the top widget
    fNewSceneTreeItemTreeWidget->insertTopLevelItem(0,item);

    // Add text that appears in the scene tree
    item->setText(0, model.GetModelType().c_str());

    // Load with info from model
    // There may be a way to add data as a QVariant, or a list of QVariants,
    // but let's try adding a G4SceneTreeItem pointer as a hex string. (There
    // does not seem to be a way of adding a pointer directly.)
    std::ostringstream oss; oss << std::hex << &model;
    auto data = QVariant(oss.str().c_str());
    item->setData(0, Qt::UserRole, data);

    // Load a tooltip
    item->setToolTip(0, model.GetModelDescription().c_str());

    // Set the check state
    item->setCheckState
    (0, model.GetVisAttributes().IsVisible()? Qt::Checked: Qt::Unchecked);

    // Set the expand state
    item->setExpanded(model.IsExpanded());

    if (model.GetType() == G4SceneTreeItem::pvmodel) {
      BuildPVQTree(model,item);
    }
  }
}

// Build Physical Volume tree of touchables
void G4UIQt::BuildPVQTree(const G4SceneTreeItem& g4stItem, QTreeWidgetItem* qtwItem)
{
  const auto& g4stChildren = g4stItem.GetChildren();
  for (const auto& g4stChild: g4stChildren) {
    QStringList qStringList;
    qStringList.append(g4stChild.GetDescription().c_str());
    auto newQTWItem = new QTreeWidgetItem(qStringList);

    // Add a GUI-side representation of the touchable as a child
    qtwItem->addChild(newQTWItem);
    
    // Load with info from g4stChild
    // There may be a way to add data as a QVariant, or a list of QVariants,
    // but let's try adding a G4SceneTreeItem pointer as a hex string. (There
    // does not seem to be a way of adding a pointer directly.)
    std::ostringstream oss; oss << std::hex << &g4stChild;
    auto data = QVariant(oss.str().c_str());
    newQTWItem->setData(0, Qt::UserRole, data);

    // Load a tooltip
    if (g4stChild.GetType() == G4SceneTreeItem::ghost) {
      auto& nameCopyNo = g4stChild.GetDescription();
      auto name = nameCopyNo.substr(0,nameCopyNo.find(':'));
      oss.str(""); oss << nameCopyNo <<
      ": Click to make visible and get more information."
      "\n  This may not work if the volume is in the \"base path\". (Hover on"
      "\n  the model to see base path.) If this is the case,"
      "\n  \"/vis/scene/add/volume " << name << "\" to bring into the displayed tree.)";
      newQTWItem->setToolTip(0, oss.str().c_str());
    } else {  // A fully defined touchable
      oss.str(""); oss << g4stChild.GetPVPath() <<
      "\nTo see properties, right-click/dump.";
      newQTWItem->setToolTip(0, oss.str().c_str());
    }

    // Set the check state
    newQTWItem->setCheckState
    (0, g4stChild.GetVisAttributes().IsVisible()? Qt::Checked: Qt::Unchecked);
    
    // Set the expand state
    newQTWItem->setExpanded(g4stChild.IsExpanded());

    // Set colour icon
    QPixmap pixmap = QPixmap(QSize(16, 16));
    pixmap.fill(ConvertG4ColourToQColor(g4stChild.GetVisAttributes().GetColour()));
    QPainter painter(&pixmap);
    painter.setPen(Qt::black);
    painter.drawRect(0,0,15,15); // Draw contour
    newQTWItem->setIcon(0,pixmap);

    // Continue recursively
    BuildPVQTree(g4stChild,newQTWItem);
  }
}

void G4UIQt::SceneTreeItemClicked(QTreeWidgetItem* item)
{
  if (item == nullptr) return;
  auto sceneTreeItem = ConvertToG4SceneTreeItem(item);
  if (sceneTreeItem == nullptr) return;

  auto uiMan = G4UImanager::GetUIpointer();

  // Respond according to type
  G4String argument = "false", inverse = "true";
  auto newCheckState = item->checkState(0);
  auto oldCheckState
  = sceneTreeItem->GetVisAttributes().IsVisible()? Qt::Checked: Qt::Unchecked;
  switch (sceneTreeItem->GetType()) {
    case G4SceneTreeItem::unidentified:
      break;  // Do nothing
    case G4SceneTreeItem::root:
      break;  // Do nothing
    case G4SceneTreeItem::model:
      // Clicked - but has checkbox actually been clicked?
      if (newCheckState != oldCheckState) {
        if (newCheckState == Qt::Checked) argument = "true";
        G4String modelName, text;
        std::istringstream iss(sceneTreeItem->GetModelDescription());
        iss >> modelName >> text;
        if (modelName.find("Text") != std::string::npos) {
          // Text model: special case, use text to identify
          uiMan->ApplyCommand("/vis/scene/activateModel " + text + ' ' + argument);
        } else {
          uiMan->ApplyCommand("/vis/scene/activateModel " + modelName + ' ' + argument);
        }
      }
      break;
    case G4SceneTreeItem::pvmodel:
      // Clicked - but has checkbox actually been clicked?
      if (newCheckState != oldCheckState) {
        if (newCheckState == Qt::Checked) argument = "true";
        G4String modelName, pvName;
        std::istringstream iss(sceneTreeItem->GetModelDescription());
        iss >> modelName >> pvName;
        uiMan->ApplyCommand("/vis/scene/activateModel " + pvName + ' ' + argument);
      }
      break;
    case G4SceneTreeItem::ghost:
      [[fallthrough]];
    case G4SceneTreeItem::touchable:
      // Construct and apply touchable commands
      // Clicked - but has checkbox actually been clicked?
      if (newCheckState != oldCheckState) {
        if (newCheckState == Qt::Checked) {
          argument = "true"; inverse = "false";
        }
        uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
        uiMan->ApplyCommand("/vis/touchable/set/visibility " + argument);
        // If daughters, set daughtersInvisible too
        if (sceneTreeItem->GetChildren().size() > 0 ) {
          uiMan->ApplyCommand("/vis/touchable/set/daughtersInvisible " + inverse);
        }
        // If not cancelled and if daughters > 0 and if making invisible
        static G4bool wanted = true;
        if (wanted && sceneTreeItem->GetChildren().size() > 0 && argument == "false") {
          QMessageBox msgBox;
          msgBox.setText
          ("This action makes this volume and all descendants invisible."
           " To see descendants, right-click and select daughtersInvisible/false"
           " and check visibility of descendants individually.");
          msgBox.setInformativeText
          ("To suppress this message click \"Discard\" or \"Don't Save\"");
          msgBox.setStandardButtons(QMessageBox::Discard | QMessageBox::Ok);
          msgBox.setDefaultButton(QMessageBox::Ok);
          auto action = msgBox.exec();
          if (action == QMessageBox::Discard) {
            wanted = false;
          }
        }
      }
      break;
  }
}

void G4UIQt::SceneTreeItemDoubleClicked(QTreeWidgetItem* item)
{
  if (item == nullptr) return;
  auto sceneTreeItem = ConvertToG4SceneTreeItem(item);
  if (sceneTreeItem == nullptr) return;
  if (sceneTreeItem->GetType() != G4SceneTreeItem::touchable) return;

  auto oldQColor = ConvertG4ColourToQColor(sceneTreeItem->GetVisAttributes().GetColour());
  auto newQColor = QColorDialog::getColor(oldQColor, fNewSceneTreeItemTreeWidget, "", QColorDialog::ShowAlphaChannel);
  if (!newQColor.isValid()) return;
  if (newQColor == oldQColor) return;

  auto newColour = ConvertQColorToG4Colour(newQColor);
  std::ostringstream oss; oss << std::setprecision(2)
  << newColour.GetRed() << ' ' << newColour.GetGreen()
  << ' ' << newColour.GetBlue() << ' ' << newColour.GetAlpha();
  auto uiMan = G4UImanager::GetUIpointer();
  uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
  uiMan->ApplyCommand("/vis/touchable/set/colour " + oss.str());
}

void G4UIQt::SceneTreeItemExpanded(QTreeWidgetItem* item)
{
  if (item == nullptr) return;
  auto sceneTreeItem = ConvertToG4SceneTreeItem(item);
  if (sceneTreeItem == nullptr) return;

  if (sceneTreeItem->GetType() == G4SceneTreeItem::ghost ||
      sceneTreeItem->GetType() == G4SceneTreeItem::touchable) {
    sceneTreeItem->SetExpanded(true);
  }
}

void G4UIQt::SceneTreeItemCollapsed(QTreeWidgetItem* item)
{
  if (item == nullptr) return;
  auto sceneTreeItem = ConvertToG4SceneTreeItem(item);
  if (sceneTreeItem == nullptr) return;

  if (sceneTreeItem->GetType() == G4SceneTreeItem::ghost ||
      sceneTreeItem->GetType() == G4SceneTreeItem::touchable) {
    sceneTreeItem->SetExpanded(false);
  }
}

void G4UIQt::NewSceneTreeItemTreeWidget::mousePressEvent(QMouseEvent* ev)
{
#if (QT_VERSION < QT_VERSION_CHECK(6, 0, 0))
  auto currentMousePressPosition = ev->globalPos();
#else
  auto currentMousePressPosition = ev->globalPosition().toPoint();
#endif

  if (ev->button() == Qt::RightButton) {
    auto item = currentItem();
    if (item) {
      auto sceneTreeItem = ConvertToG4SceneTreeItem(item);
      if (sceneTreeItem) {
        if (sceneTreeItem->GetType() == G4SceneTreeItem::touchable) {

          // We wish to present actions (and menus) in alphabetical order. I have
          // not found the required insert methods so we make our own list (map)
          // and add actions and menus and then creating appropriate qActions.
          // Using std::map to get alphabetical order. To control the order use
          // std::list of pairs instead.
          enum CommandType {withoutParameter, withABool, withAnInteger, withADouble, withAString};
          static std::map<G4String,CommandType> alphabetical;  // Ensures alphabetical order
          static G4bool first = true;
          if (first) {
            first = false;
            // Rather than simply take all commands, select those that make sense
            // in this pop-up menu.
            // Select from /vis/touchable. These are actions without parameter or for
            // which it only makes sense to use their default omitable paramater.
            alphabetical["centreAndZoomInOn"] = withoutParameter;
            alphabetical["centreOn"] = withoutParameter;
            alphabetical["dump"] = withoutParameter;
            alphabetical["extentForField"] = withoutParameter;
            alphabetical["localAxes"] = withoutParameter;
            alphabetical["showExtent"] = withoutParameter;
            alphabetical["twinkle"] = withoutParameter;
            alphabetical["volumeForField"] = withoutParameter;
            // Commands from /vis/touchable/set
            // Actions with a Boolean paramater
            alphabetical["daughtersInvisible"] = withABool;
            alphabetical["forceAuxEdgeVisible"] = withABool;
            alphabetical["forceCloud"] = withABool;
            alphabetical["forceSolid"] = withABool;
            alphabetical["forceWireframe"] = withABool;
            alphabetical["visibility"] = withABool;
            // Actions with an integer
            alphabetical["lineSegmentsPerCircle"] = withAnInteger;
            alphabetical["numberOfCloudPoints"] = withAnInteger;
            // Actions with a double
            alphabetical["lineWidth"] = withADouble;
            // Actions with a string
            alphabetical["lineStyle"] = withAString;
          }

          QMenu topMenu;  // Local (temporary) object for this item
          std::vector<QAction*> actions;  // Temporary container for action pointers
          std::vector<QMenu*> menus;  // Temporary container for menu pointers

          for (const auto& action : alphabetical) {
            const auto& af = action.first;  // G4String name of action

            if (action.second == withoutParameter) {

              auto qAction = new QAction(af.c_str(), this);
              actions.push_back(qAction);  // into temporary container
              topMenu.addAction(qAction);
              connect(qAction, &QAction::triggered, this,
                      [this, &af, &sceneTreeItem]{ActWithoutParameter(af, sceneTreeItem);});

            } else if (action.second == withABool) {

              auto menu = new QMenu(af.c_str());
              menus.push_back(menu);  // into temporary container
              topMenu.addMenu(menu);
              auto qActionTrue = new QAction("true", this);
              actions.push_back(qActionTrue);  // into temporary container
              menu->addAction(qActionTrue);
              connect(qActionTrue, &QAction::triggered, this,
                      [this, &af, &sceneTreeItem]{ActWithABool(af, sceneTreeItem, true);});
              auto qActionFalse = new QAction("false", this);
              actions.push_back(qActionFalse);  // into temporary container
              menu->addAction(qActionFalse);
              connect(qActionFalse, &QAction::triggered, this,
                      [this, &af, &sceneTreeItem]{ActWithABool(af, sceneTreeItem, false);});

            } else if (action.second == withAnInteger) {

              auto qAction = new QAction(af.c_str(), this);
              actions.push_back(qAction);  // into temporary container
              topMenu.addAction(qAction);
              connect(qAction, &QAction::triggered, this,
                      [this, &af, &sceneTreeItem]{ActWithAnInteger(af, sceneTreeItem);});

            } else if (action.second == withADouble) {

              auto qAction = new QAction(af.c_str(), this);
              actions.push_back(qAction);  // into temporary container
              topMenu.addAction(qAction);
              connect(qAction, &QAction::triggered, this,
                      [this, &af, &sceneTreeItem]{ActWithADouble(af, sceneTreeItem);});

            } else if (action.second == withAString) {

              auto qAction = new QAction(af.c_str(), this);
              actions.push_back(qAction);  // into temporary container
              topMenu.addAction(qAction);
              connect(qAction, &QAction::triggered, this,
                      [this, &af, &sceneTreeItem]{ActWithAString(af, sceneTreeItem);});
            }
          }

          topMenu.exec(currentMousePressPosition);

          // Clean up
          for (auto action : actions) {
            // No need to disconnect. Qt say, "A signal-slot connection is removed
            // when either of the objects involved are destroyed."
            delete action;
          }
          for (auto menu : menus) {
            delete menu;
          }
        }
      }
    }
  }
  
  // Pass event on up the widget tree for other actions
  QTreeWidget::mousePressEvent(ev);
}

void G4UIQt::NewSceneTreeItemTreeWidget::ActWithoutParameter
 (const G4String& action, G4SceneTreeItem* sceneTreeItem)
{
  // Special case: dump
  if (action == "dump") {
    static G4bool wanted = true;
    if (wanted) {
      QMessageBox msgBox;
      std::ostringstream oss;
      oss << G4AttCheck(sceneTreeItem->GetAttValues(), sceneTreeItem->GetAttDefs());
      // Just the first 1000 characters, otherwise it spreads off screen
      msgBox.setText((oss.str().substr(0,1000)+"...").c_str());
      msgBox.setInformativeText
      ("To suppress this message click \"Discard\" or \"Don't Save\"."
       "\nTo get a complete dump to session output click \"Ok\","
       "\nElse click \"Close\".");
      msgBox.setStandardButtons
      (QMessageBox::Discard | QMessageBox::Close | QMessageBox::Ok);
      msgBox.setDefaultButton(QMessageBox::Ok);
      auto result = msgBox.exec();
      if (result == QMessageBox::Discard) {
        wanted = false;
      } else if (result == QMessageBox::Close) {
        return;
      }
    }
  }
  auto uiMan = G4UImanager::GetUIpointer();
  uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
  uiMan->ApplyCommand("/vis/touchable/" + action);
}

void G4UIQt::NewSceneTreeItemTreeWidget::ActWithABool
 (const G4String& action, G4SceneTreeItem* sceneTreeItem, G4bool whatever)
{
  auto uiMan = G4UImanager::GetUIpointer();
  uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
  G4String which = whatever? "true": "false";
  uiMan->ApplyCommand("/vis/touchable/set/" + action + ' ' + which);
}

void G4UIQt::NewSceneTreeItemTreeWidget::ActWithAnInteger
(const G4String& action, G4SceneTreeItem* sceneTreeItem)
{
  G4bool ok = true;
  auto newValue = QInputDialog::getInt(this, action.c_str(), action.c_str(),
                                       0, 0, 999, 1, &ok);
  if (ok) {
    auto uiMan = G4UImanager::GetUIpointer();
    uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
    uiMan->ApplyCommand("/vis/touchable/set/" + action + ' '
                        + G4UIcommand::ConvertToString(newValue));
  }
}

void G4UIQt::NewSceneTreeItemTreeWidget::ActWithADouble
 (const G4String& action, G4SceneTreeItem* sceneTreeItem)
{
  G4bool ok = true;
  auto newValue = QInputDialog::getDouble(this, action.c_str(), action.c_str(),
                                          0, 0, 999, 1, &ok);
  if (ok) {
    auto uiMan = G4UImanager::GetUIpointer();
    uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
    uiMan->ApplyCommand("/vis/touchable/set/" + action + ' '
                        + G4UIcommand::ConvertToString(newValue));
  }
}

void G4UIQt::NewSceneTreeItemTreeWidget::ActWithAString
 (const G4String& action, G4SceneTreeItem* sceneTreeItem)
{
  auto uiMan = G4UImanager::GetUIpointer();
  auto command = uiMan->FindCommand("/vis/touchable/set/" + action);
  if (command) {
    QStringList qStringList;
    const auto& candidates = command->GetParameter(0)->GetParameterCandidates();
    std::istringstream iss(candidates);
    G4String candidate;
    while (iss >> candidate) qStringList.append(candidate.c_str());
    G4bool ok = true;
    auto chosenValue = QInputDialog::getItem(this, action.c_str(), action.c_str(), qStringList,
                                             0, false, &ok);
    if (ok) {
      uiMan->ApplyCommand("/vis/set/touchable" + sceneTreeItem->GetPVPath());
      G4String g4ChosenValue = chosenValue.toStdString();
      uiMan->ApplyCommand("/vis/touchable/set/" + action + ' ' + g4ChosenValue);
    }
  }
}


/** Create the History ToolBox Widget
 */
QWidget* G4UIQt::CreateHistoryTBWidget()
{
  fHistoryTBWidget = new QWidget();

  auto layoutHistoryTB = new QVBoxLayout();
  fHistoryTBTableList = new QListWidget();
  fHistoryTBTableList->setSelectionMode(QAbstractItemView::SingleSelection);
  connect(fHistoryTBTableList, SIGNAL(itemSelectionChanged()), SLOT(CommandHistoryCallback()));

  layoutHistoryTB->addWidget(fHistoryTBTableList);

  fHistoryTBWidget->setLayout(layoutHistoryTB);
  return fHistoryTBWidget;
}

/** Create the Help ToolBox Widget
 */
QWidget* G4UIQt::CreateHelpTBWidget()
{
  fHelpTBWidget = new QWidget();

  auto helpWidget = new QWidget();
  auto helpLayout = new QHBoxLayout();
  auto vLayout = new QVBoxLayout();
  fHelpVSplitter = new QSplitter(Qt::Vertical);
  fHelpLine = new QLineEdit();
  helpLayout->addWidget(new QLabel("Search :"));
  helpLayout->addWidget(fHelpLine);
  connect(fHelpLine, SIGNAL(editingFinished()), this, SLOT(LookForHelpStringCallback()));

  // Create Help tree
  FillHelpTree();

  fParameterHelpLabel = new QTextEdit();
  fParameterHelpLabel->setReadOnly(true);
  fParameterHelpTable = new QTableWidget();

  // Set layouts

  if (fHelpTreeWidget != nullptr) {
    fHelpVSplitter->addWidget(fHelpTreeWidget);
    fHelpVSplitter->setStretchFactor(0,4);
  }
  fHelpVSplitter->addWidget(fParameterHelpLabel);
  fHelpVSplitter->addWidget(fParameterHelpTable);

  fParameterHelpLabel->setVisible(false);
  fParameterHelpTable->setVisible(false);
  QSizePolicy policy = QSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);
  policy.setVerticalStretch(1);
  fParameterHelpLabel->setSizePolicy(policy);
  fParameterHelpTable->setSizePolicy(policy);

  vLayout->addWidget(helpWidget);
  vLayout->addWidget(fHelpVSplitter, 1);
  vLayout->setContentsMargins(5, 5, 5, 5);

  helpWidget->setLayout(helpLayout);
  fHelpTBWidget->setLayout(vLayout);

  return fHelpTBWidget;
}

/** Create the Cout ToolBox Widget
 */
G4UIDockWidget* G4UIQt::CreateCoutTBWidget()
{
  auto coutTBWidget = new QWidget();

  auto layoutCoutTB = new QVBoxLayout();

  fCoutTBTextArea = new QTextEdit();

  fCoutFilter = new QLineEdit();
  fCoutFilter->setToolTip("Filter output by...");

  fCoutFilter->addAction(*fSearchIcon, QLineEdit::TrailingPosition);
  fCoutFilter->setStyleSheet("border-radius:7px;");

  auto coutTBClearButton = new QPushButton();
  coutTBClearButton->setIcon(*fClearIcon);
  coutTBClearButton->setToolTip("Clear console output");
  coutTBClearButton->setStyleSheet("border-radius:7px;");
  connect(coutTBClearButton, SIGNAL(clicked()), SLOT(ClearButtonCallback()));
  connect(
    fCoutFilter, SIGNAL(textEdited(const QString&)), SLOT(CoutFilterCallback(const QString&)));

  auto coutTBSaveOutputButton = new QPushButton();
  coutTBSaveOutputButton->setIcon(*fSaveIcon);
  coutTBSaveOutputButton->setToolTip("Save console output");
  coutTBSaveOutputButton->setStyleSheet("border-radius:7px;");
  connect(coutTBSaveOutputButton, SIGNAL(clicked()), SLOT(SaveOutputCallback()));

  fCoutTBTextArea->setReadOnly(true);

  auto coutButtonWidget = new QWidget();
  auto layoutCoutTBButtons = new QHBoxLayout();

#ifdef G4MULTITHREADED
  // add all candidates to widget
  fThreadsFilterComboBox = new QComboBox();
  fThreadsFilterComboBox->setInsertPolicy(QComboBox::InsertAlphabetically);
  connect(
    fThreadsFilterComboBox, SIGNAL(activated(int)), this, SLOT(ThreadComboBoxCallback(int)));

  UpdateCoutThreadFilter();

  fThreadsFilterComboBox->setToolTip("Thread selection in output");
  layoutCoutTBButtons->addWidget(new QLabel(" Threads:"));
  layoutCoutTBButtons->addWidget(fThreadsFilterComboBox);
#endif

  layoutCoutTBButtons->addWidget(fCoutFilter);
  layoutCoutTBButtons->addWidget(coutTBClearButton);
  layoutCoutTBButtons->addWidget(coutTBSaveOutputButton);
  coutButtonWidget->setLayout(layoutCoutTBButtons);

  // reduce margins
  layoutCoutTBButtons->setContentsMargins(3, 3, 3, 0);

  layoutCoutTB->addWidget(coutButtonWidget);
  layoutCoutTB->addWidget(fCoutTBTextArea);

  coutTBWidget->setLayout(layoutCoutTB);

  fCoutTBTextArea->setMinimumSize(100, 100);

  // Command line :
  auto commandLineWidget = new QWidget();
  auto layoutCommandLine = new QHBoxLayout();

  // fill them

  fCommandLabel = new QLabel("");
  fCommandArea = new QLineEdit();

  // The QCompleter will be append at SessionStart()

  fCommandArea->activateWindow();

  fCommandArea->setFocusPolicy(Qt::StrongFocus);
  fCommandArea->setFocus(Qt::TabFocusReason);
  fCommandArea->setToolTip("Apply command");

  layoutCommandLine->addWidget(fCommandLabel);
  layoutCommandLine->addWidget(fCommandArea);

  // Connect signal
  connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));
  connect(
    fCommandArea, SIGNAL(textEdited(const QString&)), SLOT(CommandEditedCallback(const QString&)));

  commandLineWidget->setLayout(layoutCommandLine);
  commandLineWidget->setMinimumSize(50, 50);

  layoutCoutTB->addWidget(commandLineWidget);

  fCoutDockWidget = new G4UIDockWidget("Output");
  fCoutDockWidget->setAllowedAreas(Qt::TopDockWidgetArea | Qt::BottomDockWidgetArea);

  fCoutDockWidget->setWidget(coutTBWidget);
  return fCoutDockWidget;
}

/** Create the VisParameters ToolBox Widget
 */
QWidget* G4UIQt::CreateVisParametersTBWidget() { return nullptr; }

/** Create the VisParameters ToolBox Widget
 */
G4UIDockWidget* G4UIQt::CreateUITabWidget()
{
  fUITabWidget = new QTabWidget();

  // the left dock
  fUITabWidget->addTab(CreateSceneTreeWidget(), "Scene tree");
  fUITabWidget->addTab(CreateHelpTBWidget(), "Help");
  fUITabWidget->addTab(CreateHistoryTBWidget(), "History");
  fUITabWidget->setCurrentWidget(fHelpTBWidget);

  fUITabWidget->setTabToolTip(0, "Tree of scene items");
  fUITabWidget->setTabToolTip(1, "Help widget");
  fUITabWidget->setTabToolTip(2, "All commands history");
  connect(fUITabWidget, SIGNAL(currentChanged(int)), SLOT(ToolBoxActivated(int)));

  fUIDockWidget = new G4UIDockWidget("");
  fUIDockWidget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

  fUIDockWidget->setWidget(fUITabWidget);

  return fUIDockWidget;
}

QWidget* G4UIQt::CreateSceneTreeWidget()
{
  fSceneTreeWidget = new QWidget();
  auto layout = new QVBoxLayout();
  fSceneTreeWidget->setLayout(layout);

  fSceneTreeWidget->setVisible(false);

  return fSceneTreeWidget;
}

void G4UIQt::CreateViewerWidget()
{
  // Set layouts

  // clang-format off
  SetStartPage(std::string("<table width='100%'><tr><td width='30%'></td><td><div ")+
                             "style='color: rgb(140, 31, 31); font-size: xx-large; font-family: Garamond, serif; padding-bottom: 0px; font-weight: normal'>Geant4: "+
                             QApplication::applicationName ().toStdString()+
                             "</div></td><td width='40%'>&nbsp;<br/><i>http://cern.ch/geant4/</i></td></tr></table>"+
                             "<p>&nbsp;</p>"+
                             "<div style='background:#EEEEEE;'><b>Tooltips :</b><ul>"+
                             "<li><b>Start a new viewer :</b><br />"+
                             "<i>'/vis/open/...'<br />"+
                             "For example '/vis/open OGL'</i></li>"+
                             "<li><b>Execute a macro file :</b><br />"+
                             "<i>'/control/execute my_macro_file'</i></li>"+
                             "</ul></div>"+

                             "<div style='background:#EEEEEE;'><b>Documentation :</b><ul>"+
                             "<li><b>Visualisation publication :</b><br />"+
                             "<i><a href='http://www.worldscientific.com/doi/abs/10.1142/S1793962313400011'>The Geant4 Visualization System - A Multi-Driver Graphics System</b><br />,  Allison, J. et al., International Journal of Modeling, Simulation, and Scientific Computing, Vol. 4, Suppl. 1 (2013) 1340001</a>:<br/> http://www.worldscientific.com/doi/abs/10.1142/S1793962313400011</i></li>"+
                             "</ul></div>"+

                             "<div style='background:#EEEEEE;'><b>Getting Help :</b><ul>"+
                             "<li><b>If problems arise, try <a href='https://cern.ch/geant4-forum'>browsing the user forum</a> to see whether or not your problem has already been encountered.<br /> If it hasn't, you can post it and Geant4 developers will do their best to find a solution. This is also a good place to<br /> discuss Geant4 topics in general.</b> https://cern.ch/geant4-forum"+
                             "<li><b>Get a look at <a href='http://cern.ch/geant4/support'>Geant4 User support pages</a>: <i>http://cern.ch/geant4/support</i></b></li>"+
                             "</ul></div>"
                             );
  // clang-format on

  // fill right splitter
  if (fViewerTabWidget == nullptr) {
    fViewerTabWidget = new G4QTabWidget();
    fMainWindow->setCentralWidget(fViewerTabWidget);
    fViewerTabWidget->setTabsClosable(true);

    fViewerTabWidget->setUsesScrollButtons(true);

    connect(fViewerTabWidget, SIGNAL(tabCloseRequested(int)),this, SLOT(TabCloseCallback(int)));
    connect(fViewerTabWidget, SIGNAL(currentChanged(int)), SLOT(UpdateTabWidget(int)));
  }

  // set the QGLWidget size policy
  QSizePolicy policy = QSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
  policy.setVerticalStretch(4);
  fViewerTabWidget->setSizePolicy(policy);

  fViewerTabWidget->setMinimumSize(40, 40);
}

/** Get the ViewerComponents ToolBox Widget
 */
QWidget* G4UIQt::GetSceneTreeWidget() { return fSceneTreeWidget; }

/** Get the Viewer properties  Widget
 */
QWidget* G4UIQt::GetViewerPropertiesWidget()
{
  if (fViewerPropertiesDialog == nullptr) {
    CreateViewerPropertiesDialog();
  }
  return fViewerPropertiesWidget;
}

/** Get the Pick Widget
 */
QWidget* G4UIQt::GetPickInfosWidget()
{
  if (fPickInfosDialog == nullptr) {
    CreatePickInfosDialog();
  }
  return fPickInfosWidget;
}

/**   Add a new tab in the viewer
 */
G4bool G4UIQt::AddViewerTab(QWidget* aWidget, std::string title)
{
  if (fViewerTabWidget == nullptr) {
    return false;
  }
  fViewerTabWidget->addTab(aWidget, title.c_str());

  return true;
}

/**   Add a new tab in the viewer
 */
G4bool G4UIQt::AddViewerTabFromFile(std::string fileName, std::string title)
{
  if (fViewerTabWidget == nullptr) {
    return false;
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return false;
  std::ifstream file(UI->FindMacroPath(fileName.c_str()).data());
  if (file) {
    std::string content((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

    auto text = new QTextEdit();
    text->setAcceptRichText(true);
    text->setContentsMargins(5, 5, 5, 5);
    text->setText(QString("<pre>") + content.c_str() + "</pre>");
    text->setReadOnly(true);
    fViewerTabWidget->addTab(text, title.c_str());
  }
  else {
    return false;
  }
  return true;
}

/**   Add a new tab widget.
  Create the tab if it was not done
*/
G4bool G4UIQt::AddTabWidget(QWidget* aWidget, QString name)
{
  if (fViewerTabWidget == nullptr) {
    CreateViewerWidget();
  }

  if (aWidget == nullptr) {
    return false;
  }
  // Has to be added before we put it into the fViewerTabWidget widget
  aWidget->setParent(fViewerTabWidget);  // Will create in some cases widget outside
  // of UI for a really short moment

  fViewerTabWidget->addTab(aWidget, name);

  fViewerTabWidget->setCurrentIndex(fViewerTabWidget->count() - 1);

  // Set visible
  fViewerTabWidget->setLastTabCreated(fViewerTabWidget->currentIndex());

  // Not the good solution, but ensure that the help tree is correctly build when launching a viewer
  // It should be done by a notification when adding a command, but that's nit done yet
  // (Geant4.10.1)
  FillHelpTree();

  return true;
}

void G4UIQt::SetStartPage(const std::string& text)
{
  if (! text.empty()) {
    fDefaultViewerFirstPageHTMLText = text;
  }
  if (fStartPage == nullptr) {
    fStartPage = new QTextBrowser();
    fStartPage->setContentsMargins(5, 5, 5, 5);
    fStartPage->setReadOnly(true);
  }
  fStartPage->setOpenExternalLinks(true);
  fStartPage->setHtml(fDefaultViewerFirstPageHTMLText.c_str());
}

void G4UIQt::UpdateTabWidget(int tabNumber)
{
  if (fViewerTabWidget == nullptr) {
    fViewerTabWidget = new G4QTabWidget;
  }

  fViewerTabWidget->setCurrentIndex(tabNumber);

  // Send this signal to unblock graphic updates !
  fViewerTabWidget->setTabSelected(false);

  fViewerTabWidget->setVisible(true);

  // This will send a paintEvent to OGL Viewers
  fViewerTabWidget->setTabSelected(true);
}

/** Send resize event to all tabs
 */
void G4UIQt::ResizeTabWidget(QResizeEvent* e)
{
  if (fViewerTabWidget != nullptr) {
    for (G4int a = 0; a < fViewerTabWidget->count(); a++) {
      fViewerTabWidget->widget(a)->resize(e->size());
    }
  }
}

/**   Start the Qt main loop
 */
G4UIsession* G4UIQt::SessionStart()
{
  G4Qt* interactorManager = G4Qt::getInstance();
  Prompt("Session :");
  exitSession = false;

  QCoreApplication::sendPostedEvents();

  fMainWindow->setVisible(true);

  if (fDefaultIcons) {
    fToolbarApp->setVisible(true);
  }
  else {
    // Set not visible until session start
    fToolbarApp->setVisible(false);
  }
  // Rebuild help tree (new command could be registered)
  FillHelpTree();

  // Rebuild command completion (new command could be registered)
  UpdateCommandCompleter();

  // Set event filters
  fHistoryTBTableList->installEventFilter(this);
  fCommandArea->installEventFilter(this);

  // Focus on command line
  fCommandArea->setFocus();

  interactorManager->DisableSecondaryLoop();  // TO KEEP
  if ((QApplication*)interactorManager->GetMainInteractor() != nullptr)
    ((QApplication*)interactorManager->GetMainInteractor())->exec();

  interactorManager->EnableSecondaryLoop();
  return this;
}

/**   Display the prompt in the prompt area
   @param aPrompt : string to display as the promt label
*/
void G4UIQt::Prompt(G4String aPrompt)
{
  if (aPrompt == nullptr) return;

  fCommandLabel->setText((char*)aPrompt.data());
}

void G4UIQt::SessionTerminate()
{
  G4Qt* interactorManager = G4Qt::getInstance();
  fMainWindow->close();
  ((QApplication*)interactorManager->GetMainInteractor())->exit();
}

/**
   Called by intercoms/src/G4UImanager.cc<br>
   Called by visualization/management/src/G4VisCommands.cc with "EndOfEvent" argument<br>
   It have to pause the session command terminal.<br>
   Call SecondaryLoop to wait for exit event<br>
   @param aState
   @see : G4VisCommandReviewKeptEvents::SetNewValue
*/
void G4UIQt::PauseSessionStart(const G4String& aState)
{
  if (aState == nullptr) return;

  if (aState == "G4_pause> ") {  // TO KEEP
    SecondaryLoop("Pause, type continue to exit this state");  // TO KEEP
  }  // TO KEEP

  if (aState == "EndOfEvent") {  // TO KEEP
    // Picking with feed back in event data Done here !!!
    SecondaryLoop("End of event, type continue to exit this state");  // TO KEEP
  }  // TO KEEP
}

/**
   Begin the secondary loop
   @param a_prompt : label to display as the prompt label
 */
void G4UIQt::SecondaryLoop(G4String aPrompt)
{
  if (aPrompt == nullptr) return;

  G4Qt* interactorManager = G4Qt::getInstance();  // TO KEEP ?
  Prompt(aPrompt);  // TO KEEP
  exitPause = false;  // TO KEEP
  while (true) {
    ((QApplication*)interactorManager)->processEvents(QEventLoop::WaitForMoreEvents);
    if (exitPause) break;  // TO KEEP
  }  // TO KEEP
  Prompt("Session :");  // TO KEEP
}

#ifdef G4MULTITHREADED
#  include "G4AutoLock.hh"
#  include "G4Threading.hh"
namespace
{
  G4Mutex ReceiveMutex = G4MUTEX_INITIALIZER;
}  // namespace
#endif

/**
   Receive a debug log message from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIQt::ReceiveG4debug(const G4String& aString)
{
  if (aString.empty()) return 0;

#ifdef G4MULTITHREADED
  G4AutoLock al(&ReceiveMutex);
#endif

  // A workaround so that output is not lost after crash or G4Exception.
  // The "workaround" is to make sure all flushed output appears on
  // the terminal after a crash, because even flushed output can
  // get lost in the Qt UI system.
  // But...it seems workers write to std::cout/cerr anyway (is that a bug?),
  // so limit this to the master thread
#ifdef G4MULTITHREADED
  if (G4Threading::IsMasterThread())
#endif
    std::cout << aString << std::flush;

  G4String aStringWithStyle;
  // aString has a \n on the end (maybe it comes from G4endl or from the
  // Enter key on the command line) - ignore it. That’s why
  // i < aString.length() - 1
  // But other \n need to be translated to an HTML newline.
  // Similarly, spaces need to be translated to an HTML "non-breaking space".
  // Tabs (\t) are more tricky since the number of equivalent spaces depends
  // on how many characters precede it. Probably needs an HTML table. For now
  // we replace \t with four spaces.
  for (G4int i = 0; i < (G4int)aString.length() - 1; ++i) {
    if (aString[i] == '\n') {
      aStringWithStyle += "<br>";
    }
    else if (aString[i] == ' ') {
      aStringWithStyle += "&nbsp;";
    }
    else if (aString[i] == '\t') {
      aStringWithStyle += "&nbsp;&nbsp;&nbsp;&nbsp;";
    }
    else if (aString[i] == '<') {
      aStringWithStyle += "&lt;";
    }
    else {
      aStringWithStyle += aString[i];
    }
  }
  if (fOutputStyles["debug"].fixed) {
    aStringWithStyle = "<span style='font-family:courier;'>" + aStringWithStyle + "</span>";
  }
  else {
    aStringWithStyle = "<span>" + aStringWithStyle + "</span>";
  }

  // Add to string
  G4UIOutputString txt =
    G4UIOutputString(QString((char*)aStringWithStyle.data()), GetThreadPrefix());
  fG4OutputString.push_back(txt);

#ifdef G4MULTITHREADED
  QString result = FilterOutput(txt, fThreadsFilterComboBox->currentText(), fCoutFilter->text());
#else
  QString result = FilterOutput(txt, "", fCoutFilter->text());
#endif

  if (result.isEmpty()) {
    return 0;
  }

  if (fOutputStyles["debug"].highlight) {
    QPalette pal;
    result = QString("<span style='background:") + pal.link().color().name() + ";'>&nbsp;</span>" +
             "<span style='background: Pink;'> " + result + "</span>";
  }
  result = QString("<font color=\"Green\">") + result + QString("</font>");

  fCoutTBTextArea->append(result);
  fCoutTBTextArea->ensureCursorVisible();

#ifdef G4MULTITHREADED
  UpdateCoutThreadFilter();
#endif

  return 0;
}

/**
   Receive a cout from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIQt::ReceiveG4cout(const G4String& aString)
{
  if (aString.empty()) return 0;

  // Try to be smart :
  // "*** This is just a warning message. ***"
  if (G4StrUtil::contains(aString, "*** This is just a warning message. ***")) {
    return ReceiveG4cerr(aString);
  }

#ifdef G4MULTITHREADED
  G4AutoLock al(&ReceiveMutex);
#endif

  // A workaround so that output is not lost after crash or G4Exception.
  // The "workaround" is to make sure all flushed output appears on
  // the terminal after a crash, because even flushed output can
  // get lost in the Qt UI system.
  // But...it seems workers write to std::cout/cerr anyway (is that a bug?),
  // so limit this to the master thread
#ifdef G4MULTITHREADED
  if (G4Threading::IsMasterThread())
#endif
    std::cout << aString << std::flush;

  G4String aStringWithStyle;
  // aString has a \n on the end (maybe it comes from G4endl or from the
  // Enter key on the command line) - ignore it. That’s why
  // i < aString.length() - 1
  // But other \n need to be translated to an HTML newline.
  // Similarly, spaces need to be translated to an HTML "non-breaking space".
  // Tabs (\t) are more tricky since the number of equivalent spaces depends
  // on how many characters precede it. Probably needs an HTML table. For now
  // we replace \t with four spaces.
  for (G4int i = 0; i < (G4int)aString.length() - 1; ++i) {
    if (aString[i] == '\n') {
      aStringWithStyle += "<br>";
    }
    else if (aString[i] == ' ') {
      aStringWithStyle += "&nbsp;";
    }
    else if (aString[i] == '\t') {
      aStringWithStyle += "&nbsp;&nbsp;&nbsp;&nbsp;";
    }
    else if (aString[i] == '<') {
      aStringWithStyle += "&lt;";
    }
    else {
      aStringWithStyle += aString[i];
    }
  }
  if (fOutputStyles["cout"].fixed) {
    aStringWithStyle = "<span style='font-family:courier;'>" + aStringWithStyle + "</span>";
  }
  else {
    aStringWithStyle = "<span>" + aStringWithStyle + "</span>";
  }

  // Add to string
  G4UIOutputString txt =
    G4UIOutputString(QString((char*)aStringWithStyle.data()), GetThreadPrefix());
  fG4OutputString.push_back(txt);

#ifdef G4MULTITHREADED
  QString result = FilterOutput(txt, fThreadsFilterComboBox->currentText(), fCoutFilter->text());
#else
  QString result = FilterOutput(txt, "", fCoutFilter->text());
#endif

  if (result.isEmpty()) {
    return 0;
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (fOutputStyles["cout"].highlight) {
    if (! UI->IsLastCommandOutputTreated()) {
      QPalette pal;
      result = QString("<span style='background:") + pal.link().color().name() +
               ";'>&nbsp;</span>" + "<span style='background:" + pal.highlight().color().name() +
               ";'> " + result + "</span>";
    }
  }
  UI->SetLastCommandOutputTreated();

  fCoutTBTextArea->append(result);
  fCoutTBTextArea->ensureCursorVisible();

#ifdef G4MULTITHREADED
  UpdateCoutThreadFilter();
#endif

  // reset error stack
  fLastErrMessage = aString;
  return 0;
}

/**
   Receive a cerr from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIQt::ReceiveG4cerr(const G4String& aString)
{
  if (aString.empty()) return 0;

#ifdef G4MULTITHREADED
  G4AutoLock al(&ReceiveMutex);
#endif

  // A workaround so that output is not lost after crash or G4Exception.
  // The "workaround" is to make sure all flushed output appears on
  // the terminal after a crash, because even flushed output can
  // get lost in the Qt UI system.
  // But...it seems workers write to std::cout/cerr anyway (is that a bug?),
  // so limit this to the master thread
#ifdef G4MULTITHREADED
  if (G4Threading::IsMasterThread())
#endif
    std::cerr << aString << std::flush;

  G4String aStringWithStyle;
  // aString has a \n on the end (maybe it comes from G4endl or from the
  // Enter key on the command line) - ignore it. That’s why
  // i < aString.length() - 1
  // But other \n need to be translated to an HTML newline.
  // Similarly, spaces need to be translated to an HTML "non-breaking space".
  // Tabs (\t) are more tricky since the number of equivalent spaces depends
  // on how many characters precede it. Probably needs an HTML table. For now
  // we replace \t with four spaces.
  for (G4int i = 0; i < (G4int)aString.length() - 1; ++i) {
    if (aString[i] == '\n') {
      aStringWithStyle += "<br>";
    }
    else if (aString[i] == ' ') {
      aStringWithStyle += "&nbsp;";
    }
    else if (aString[i] == '\t') {
      aStringWithStyle += "&nbsp;&nbsp;&nbsp;&nbsp;";
    }
    else if (aString[i] == '<') {
      aStringWithStyle += "&lt;";
    }
    else {
      aStringWithStyle += aString[i];
    }
  }
  if (fOutputStyles["cerr"].fixed) {
    aStringWithStyle = "<span style='font-family:courier;'>" + aStringWithStyle + "</span>";
  }
  else {
    aStringWithStyle = "<span>" + aStringWithStyle + "</span>";
  }

  // Add to string

  G4UIOutputString txt =
    G4UIOutputString(QString((char*)aStringWithStyle.data()).trimmed(), GetThreadPrefix(), "error");
  fG4OutputString.push_back(txt);

#ifdef G4MULTITHREADED
  QString result = FilterOutput(txt, fThreadsFilterComboBox->currentText(), fCoutFilter->text());
#else
  QString result = FilterOutput(txt, "", fCoutFilter->text());
#endif
  if (result.isEmpty()) {
    return 0;
  }

  // Suppress space, \n,\t,\r...
  if (QString(aString.data()).trimmed() != "") {
    if ((G4StateManager::GetStateManager()->GetCurrentState() == G4State_Abort) ||
        (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Quit))
    {
      // In case of Abort or Quit, the useful error message should be in the last error message !
      fLastErrMessage += "\n" + aString;
      QString criticalMessage = fLastErrMessage.data();
      criticalMessage = criticalMessage.toHtmlEscaped();
      QMessageBox::critical(fMainWindow, "Error", QString(fLastErrMessage));
    }
  }
  fCoutTBTextArea->append(QString("<font color=\"Red\">") + result + QString("</font>"));
  fCoutTBTextArea->ensureCursorVisible();

  if (QString(aString.data()).trimmed() != "") {
    fLastErrMessage += aString;
  }
#ifdef G4MULTITHREADED
  UpdateCoutThreadFilter();
#endif
  return 0;
}

G4String G4UIQt::GetThreadPrefix()
{
  G4String threadPrefix = "";
#ifdef G4MULTITHREADED
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return "";
  if (UI->GetThreadCout() != nullptr) {
    threadPrefix = UI->GetThreadCout()->GetFullPrefixString().data();
    if (UI->GetThreadCout()->GetPrefixString() == G4String("G4VIS")) {
      return "G4VIS";
    }
  }
#endif
  return threadPrefix;
}

#ifdef G4MULTITHREADED
void G4UIQt::UpdateCoutThreadFilter()
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return;

  // add "All" and "Master"
  if (fThreadsFilterComboBox->count() < 2) {
    if (fThreadsFilterComboBox->findText("All", Qt::MatchExactly) == -1) {
      fThreadsFilterComboBox->addItem("All");
    }
  }
  if (fThreadsFilterComboBox->count() < 2) {
    if (fThreadsFilterComboBox->findText("Master", Qt::MatchExactly) == -1) {
      fThreadsFilterComboBox->addItem("Master");
    }
  }
  // Add current Cout
  G4String prefix = GetThreadPrefix();
  if (! prefix.empty()) {
    if (fThreadsFilterComboBox->findText(prefix.data(), Qt::MatchExactly) == -1) {
      fThreadsFilterComboBox->addItem(prefix.data());
    }
  }
}
#endif

/**
   Add a new menu to the menu bar
   @param aName name of menu
   @param aLabel label to display
 */
void G4UIQt::AddMenu(const char* aName, const char* aLabel)
{
  if (aName == nullptr) return;
  if (aLabel == nullptr) return;

  auto fileMenu = new QMenu(aLabel);
  fMainWindow->menuBar()->addMenu(fileMenu);

  AddInteractor(aName, (G4Interactor)fileMenu);
}

/**
   Add a new button to a menu
   @param aMenu : parent menu
   @param aLabel : label to display
   @param aCommand : command to execute as a callback
 */
void G4UIQt::AddButton(const char* aMenu, const char* aLabel, const char* aCommand)
{
  if (aMenu == nullptr) return;  // TO KEEP
  if (aLabel == nullptr) return;  // TO KEEP
  if (aCommand == nullptr) return;  // TO KEEP

  QMenu* parentTmp = (QMenu*)GetInteractor(aMenu);

  if (parentTmp == nullptr) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();

    if (verbose >= 2) {
      G4cout << "Menu name " << aMenu << " does not exist, please define it before using it."
             << G4endl;
    }
    return;
  }

  // Find the command in the command tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return;
  G4UIcommandTree* treeTop = UI->GetTree();

  G4String cmd = aCommand;
  std::size_t cmdEndPos = cmd.find_first_of(" \t");
  if (cmdEndPos != std::string::npos) {
    cmd.erase(cmdEndPos);
  }

  if (treeTop->FindPath(cmd) == nullptr) {
    if (cmd != "ls" && cmd.substr(0, 3) != "ls " && cmd != "pwd" && cmd != "cd" &&
        cmd.substr(0, 3) != "cd " && cmd != "help" && cmd.substr(0, 5) != "help " &&
        cmd[0] != '?' && cmd != "hist" && cmd != "history" && cmd[0] != '!' && cmd != "exit" &&
        cmd != "cont" && cmd != "continue")
    {
      G4UImanager* UImanager = G4UImanager::GetUIpointer();
      G4int verbose = UImanager->GetVerboseLevel();

      if (verbose >= 2) {
        G4cout << "Warning: command '" << cmd
               << "' does not exist, please define it before using it." << G4endl;
      }
    }
  }

  QString cmd_tmp = QString(aCommand);
  parentTmp->addAction(aLabel, this, [this, cmd_tmp]() { this->ButtonCallback(cmd_tmp); });
}

/**
 special case for the "open" icon. It will open a file selector and map the return file to the given
 command.
*/
void G4UIQt::AddIcon(
  const char* aLabel, const char* aIconFile, const char* aCommand, const char* aFileName)
{
  if (aLabel == nullptr) return;  // TO KEEP
  // special case, aCommand could be NULL if aIconFile is not user_icon
  if (aCommand == nullptr) {
    if (std::string(aIconFile) == "user_icon") {
      return;  // TO KEEP
    }
  }
  QPixmap* pix;
  G4bool userToolBar = false;

  if (! fDefaultIcons) {
    userToolBar = true;
  }
  if (std::string(aIconFile) == "user_icon") {
    // try to open a file
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    pix = new QPixmap(UImanager->FindMacroPath(aFileName).data());
    if (pix->isNull()) {
      G4int verbose = UImanager->GetVerboseLevel();

      if (verbose >= 2) {
        G4cout << "Warning: file '" << aFileName
               << "' is incorrect or does not exist, this command will not be build" << G4endl;
      }
      return;
    }
  }
  else if (std::string(aIconFile) == "open") {
    pix = fOpenIcon;
  }
  else if (std::string(aIconFile) == "save") {
    pix = fSaveIcon;
  }
  else if (std::string(aIconFile) == "move") {
    pix = fMoveIcon;
  }
  else if (std::string(aIconFile) == "rotate") {
    pix = fRotateIcon;
  }
  else if (std::string(aIconFile) == "pick") {
    pix = fPickIcon;
  }
  else if (std::string(aIconFile) == "zoom_in") {
    pix = fZoomInIcon;
  }
  else if (std::string(aIconFile) == "zoom_out") {
    pix = fZoomOutIcon;
  }
  else if (std::string(aIconFile) == "wireframe") {
    pix = fWireframeIcon;
  }
  else if (std::string(aIconFile) == "solid") {
    pix = fSolidIcon;
  }
  else if (std::string(aIconFile) == "hidden_line_removal") {
    pix = fHiddenLineRemovalIcon;
  }
  else if (std::string(aIconFile) == "hidden_line_and_surface_removal") {
    pix = fHiddenLineAndSurfaceRemovalIcon;
  }
  else if (std::string(aIconFile) == "perspective") {
    pix = fPerspectiveIcon;
  }
  else if (std::string(aIconFile) == "ortho") {
    pix = fOrthoIcon;
  }
  else if (std::string(aIconFile) == "runBeamOn") {
    pix = fRunIcon;
  }
  else if (std::string(aIconFile) == "exit") {
    pix = fExitIcon;
  }
  else {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();

    if (verbose >= 2) {
      G4cout << "Parameter" << aIconFile << " not defined" << G4endl;
    }
    return;
  }
  QToolBar* currentToolbar = nullptr;
  if (userToolBar) {
    if (fToolbarUser == nullptr) {
      fToolbarUser = new QToolBar();
      fToolbarUser->setIconSize(QSize(20, 20));
      fMainWindow->addToolBar(Qt::TopToolBarArea, fToolbarUser);
    }
    currentToolbar = fToolbarUser;
  }
  else {
    if (fToolbarApp == nullptr) {
      fToolbarApp = new QToolBar();
      fToolbarApp->setIconSize(QSize(20, 20));
      fMainWindow->addToolBar(Qt::TopToolBarArea, fToolbarApp);
    }
    currentToolbar = fToolbarApp;
  }

  // Check if already present

  QList<QAction*> list = currentToolbar->actions();

  for (auto i : list) {
    if (i->text() == QString(aLabel)) {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      if (UI == nullptr) return;
      G4int verbose = UI->GetVerboseLevel();
      if (verbose >= 2) {
        G4cout << "Warning: A toolBar icon \"" << aLabel << "\" already exists with the same name!"
               << G4endl;
      }
    }
  }

  // special cases :"open"
  if (std::string(aIconFile) == "open") {
    QString txt = aCommand + fStringSeparator + aLabel;
    currentToolbar->addAction(
      QIcon(*pix), aIconFile, this, [this, txt]() { this->OpenIconCallback(txt); });

    // special cases :"save"
  }
  else if (std::string(aIconFile) == "save") {
    QString txt = aCommand + fStringSeparator + aLabel;
    currentToolbar->addAction(
      QIcon(*pix), aIconFile, this, [this, txt]() { this->SaveIconCallback(txt); });
    // special cases : cursor style
  }
  else if ((std::string(aIconFile) == "move") || (std::string(aIconFile) == "rotate") ||
           (std::string(aIconFile) == "pick") || (std::string(aIconFile) == "zoom_out") ||
           (std::string(aIconFile) == "zoom_in"))
  {
    QString txt = QString(aIconFile);
    QAction* action = currentToolbar->addAction(
      QIcon(*pix), aIconFile, this, [this, txt]() { this->ChangeCursorAction(txt); });
    action->setCheckable(true);
    action->setChecked(true);
    action->setData(aIconFile);

    if (std::string(aIconFile) == "move") {
      SetIconMoveSelected();
    }
    if (std::string(aIconFile) == "rotate") {
      SetIconRotateSelected();
    }
    if (std::string(aIconFile) == "pick") {
      SetIconPickSelected();
    }
    if (std::string(aIconFile) == "zoom_in") {
      SetIconZoomInSelected();
    }
    if (std::string(aIconFile) == "zoom_out") {
      SetIconZoomOutSelected();
    }

    // special case : surface style
  }
  else if ((std::string(aIconFile) == "hidden_line_removal") ||
           (std::string(aIconFile) == "hidden_line_and_surface_removal") ||
           (std::string(aIconFile) == "solid") || (std::string(aIconFile) == "wireframe"))
  {
    QString txt = QString(aIconFile);
    QAction* action = currentToolbar->addAction(
      QIcon(*pix), aIconFile, this, [this, txt]() { this->ChangeSurfaceStyle(txt); });
    action->setCheckable(true);
    action->setChecked(true);
    action->setData(aIconFile);

    if (std::string(aIconFile) == "hidden_line_removal") {
      SetIconHLRSelected();
    }
    if (std::string(aIconFile) == "hidden_line_and_surface_removal") {
      SetIconHLHSRSelected();
    }
    if (std::string(aIconFile) == "solid") {
      SetIconSolidSelected();
    }
    if (std::string(aIconFile) == "wireframe") {
      SetIconWireframeSelected();
    }

    // special case : perspective/ortho
  }
  else if ((std::string(aIconFile) == "perspective") || (std::string(aIconFile) == "ortho")) {
    QString txt = QString(aIconFile);
    QAction* action = currentToolbar->addAction(
      QIcon(*pix), aIconFile, this, [this, txt]() { this->ChangePerspectiveOrtho(txt); });
    action->setCheckable(true);
    action->setChecked(true);
    action->setData(aIconFile);

    if (std::string(aIconFile) == "perspective") {
      SetIconPerspectiveSelected();
    }
    if (std::string(aIconFile) == "ortho") {
      SetIconOrthoSelected();
    }
  }
  else {
    // Find the command in the command tree
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI == nullptr) return;
    G4UIcommandTree* treeTop = UI->GetTree();
    if (aCommand != nullptr) {
      std::string str = aCommand;
      std::string::size_type pos = str.find(' ');
      if (pos != std::string::npos) {
        str = str.substr(0, pos).c_str();
      }
      if (treeTop->FindPath(str.c_str()) == nullptr) {
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        G4int verbose = UImanager->GetVerboseLevel();

        if (verbose >= 2) {
          G4cout << "Warning: command '" << aCommand
                 << "' does not exist, please define it before using it." << G4endl;
        }
      }
    }
    QString txt = QString(aCommand);
    currentToolbar->addAction(
      QIcon(*pix), aCommand, this, [this, txt]() { this->ButtonCallback(txt); });
  }
}

void G4UIQt::SetOutputStyle(const char* destination, const char* style)
{
  // Specify an output style
  // First argument destination ("cout" etc or "all")
  // Second argument is the required style - see guidance

  SetStyleUtility(destination, style);
}

void G4UIQt::NativeMenu(G4bool aVal)
{
  if (fMainWindow->menuBar()->isNativeMenuBar() == aVal) return;  // already in this state

  // Menu become empty when goin from Qt to Native Bar
  fMainWindow->menuBar()->setNativeMenuBar(aVal);
}

void G4UIQt::ClearMenu() { fMainWindow->menuBar()->clear(); }

void G4UIQt::ActivateCommand(G4String newCommand)
{
  if (fHelpTreeWidget == nullptr) {
    return;
  }
  // Look for the choosen command "newCommand"
  std::size_t i = newCommand.find(' ');
  G4String targetCom = "";
  if (i != std::string::npos) {
    G4String newValue = newCommand.substr(i + 1, newCommand.length() - (i + 1));
    G4StrUtil::strip(newValue);
    targetCom = ModifyToFullPathCommand(newValue);
  }
  if (! targetCom.empty()) {
    OpenHelpTreeOnCommand(targetCom.data());
  }

  fUITabWidget->setCurrentWidget(fHelpTBWidget);
}

/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIQt::InitHelpTreeAndVisParametersWidget()
{
  if (fHelpTreeWidget == nullptr) {
    fHelpTreeWidget = new QTreeWidget();
  }

  // build widget
  fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  QStringList labels;
  labels << QString("Command");
  fHelpTreeWidget->setHeaderLabels(labels);

  connect(fHelpTreeWidget, SIGNAL(itemSelectionChanged()), this, SLOT(HelpTreeClicCallback()));
  connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), this,
    SLOT(HelpTreeDoubleClicCallback()));
}
/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIQt::FillHelpTree()
{
  if (fHelpTreeWidget == nullptr) {
    InitHelpTreeAndVisParametersWidget();
  }

  QString searchText = fHelpLine->text();

  if (searchText == "") {
    // clear old help tree
    //    fHelpTreeWidget->clear();
  }
  else {
    return;
  }

  if (fParameterHelpLabel != nullptr) {
    fParameterHelpLabel->setText("Choose a command in the command tree");
    fParameterHelpTable->setVisible(false);
  }

  if (fHelpLine != nullptr) {
    fHelpLine->setText("");
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return;
  G4UIcommandTree* treeTop = UI->GetTree();

  G4int treeSize = treeTop->GetTreeEntry();
  QTreeWidgetItem* newItem = nullptr;
  QString commandText = "";
  for (G4int a = 0; a < treeSize; ++a) {
    // Creating new item
    newItem = nullptr;

    commandText = QString((char*)(treeTop->GetTree(a + 1)->GetPathName()).data()).trimmed();

    // if already exist, don't create it !
    for (G4int b = 0; b < fHelpTreeWidget->topLevelItemCount(); ++b) {
      if (newItem == nullptr) newItem = FindTreeItem(fHelpTreeWidget->topLevelItem(b), commandText);
    }

    if (newItem == nullptr) {
      newItem = new QTreeWidgetItem();
      newItem->setText(0, GetShortCommandPath(commandText));
      fHelpTreeWidget->addTopLevelItem(newItem);
    }

    // look for childs
    CreateHelpTree(newItem, treeTop->GetTree(a + 1));
  }
}

/**   Fill the Help Tree Widget
   @param aParent : parent item to fill
   @param aCommandTree : commandTree node associate with this part of the Tree
*/
void G4UIQt::CreateHelpTree(QTreeWidgetItem* aParent, G4UIcommandTree* aCommandTree)
{
  if (aParent == nullptr) return;
  if (aCommandTree == nullptr) return;

  // Creating new item
  QTreeWidgetItem* newItem;

  QString commandText = "";
  // Get the Sub directories
  for (G4int a = 0; a < aCommandTree->GetTreeEntry(); ++a) {
    commandText = QString((char*)(aCommandTree->GetTree(a + 1)->GetPathName()).data()).trimmed();

    // if already exist, don't create it !
    newItem = FindTreeItem(aParent, commandText);
    if (newItem == nullptr) {
      newItem = new QTreeWidgetItem();
      newItem->setText(0, GetShortCommandPath(commandText));
      aParent->addChild(newItem);
    }
    CreateHelpTree(newItem, aCommandTree->GetTree(a + 1));
  }

  // Get the Commands

  for (G4int a = 0; a < aCommandTree->GetCommandEntry(); ++a) {
    QStringList stringList;
    commandText =
      QString((char*)(aCommandTree->GetCommand(a + 1)->GetCommandPath()).data()).trimmed();

    // if already exist, don't create it !
    newItem = FindTreeItem(aParent, commandText);
    if (newItem == nullptr) {
      newItem = new QTreeWidgetItem();
      newItem->setText(0, GetShortCommandPath(commandText));
      aParent->addChild(newItem);
      newItem->setExpanded(false);
    }
  }
}

/**
 Add the following command to the corresponding groupbox
 If depthLevel is 1 : create ToolBox
 If depthLevel is 2 or more : create GroupBox
*/
G4bool G4UIQt::CreateVisCommandGroupAndToolBox(
  G4UIcommand* aCommand, QWidget* aParent, G4int aDepthLevel, G4bool isDialog)
{
  QString commandText =
    QString((char*)(aCommand->GetCommandPath().data())).section("/", -aDepthLevel);

  if (commandText == nullptr) {
    return false;
  }

  // Look if groupBox is create
  //  QGroupBox* gBoxCommandWidget;
  QWidget* newParentWidget = nullptr;
  G4bool found = false;
  QString commandSection = commandText.left(commandText.indexOf("/"));

  if (aDepthLevel == 1) {
    auto currentParent = dynamic_cast<QToolBox*>(aParent);
    if (currentParent != nullptr) {
      // already exists ?
      for (G4int a = 0; a < currentParent->count(); ++a) {
        if (currentParent->itemText(a) == commandSection) {
          found = true;
          newParentWidget = currentParent->widget(a);
        }
      }
    }
    // Not found ? create it
    if (! found) {
      newParentWidget = new QGroupBox();
      newParentWidget->setLayout(new QVBoxLayout());
      if (currentParent != nullptr) {
        currentParent->addItem(newParentWidget, commandSection);
      }
      else {
        if (aParent->layout() == nullptr) {
          aParent->setLayout(new QVBoxLayout());
        }
        aParent->layout()->addWidget(newParentWidget);
      }

      if (commandText.indexOf("/") == -1) {
        // Guidance
        QString guidance;
        auto n_guidanceEntry = (G4int)aCommand->GetGuidanceEntries();
        for (G4int i_thGuidance = 0; i_thGuidance < n_guidanceEntry; i_thGuidance++) {
          guidance += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
        }
        newParentWidget->setToolTip(guidance);
      }

      auto sc = dynamic_cast<QScrollArea*>(newParentWidget->parent()->parent());
      if (sc != nullptr) {
        sc->ensureWidgetVisible(newParentWidget);
      }
    }
  }
  else {
    // try to know if this level is already there
    auto currentParent = dynamic_cast<QGroupBox*>(aParent);
    if (currentParent != nullptr) {
      // if depth==2, then we add a [more parameters inside] to the toolBoxItem parent
      // QGroupBox > QWidget > QScrollArea > QToolBox
      if (aDepthLevel == 2) {
        auto parentToolBox = dynamic_cast<QToolBox*>(currentParent->parent()->parent()->parent());
        if (parentToolBox != nullptr) {
          //          parentToolBox->setItemText(parentToolBox->indexOf(currentParent),"[more
          //          parameters inside]");
        }
      }
      for (G4int a = 0; a < aParent->layout()->count(); ++a) {
        auto gb = dynamic_cast<QGroupBox*>(aParent->layout()->itemAt(a)->widget());
        if (gb != nullptr) {
          if (gb->title() == commandSection) {
            found = true;
            newParentWidget = gb;
          }
        }
      }
    }

    // Not found ? create it
    if (! found) {
      newParentWidget = new QGroupBox();
      newParentWidget->setLayout(new QVBoxLayout());
      if (aParent->layout() == nullptr) {
        aParent->setLayout(new QVBoxLayout());
      }
      aParent->layout()->addWidget(newParentWidget);

      // set toolTip
      // Guidance
      QString guidance;
      auto n_guidanceEntry = (G4int)aCommand->GetGuidanceEntries();
      for (G4int i_thGuidance = 0; i_thGuidance < n_guidanceEntry; i_thGuidance++) {
        guidance += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
      }
      newParentWidget->setToolTip(guidance);
    }
  }

  // fill command groupbox
  if (commandText.indexOf("/") == -1) {
    if (CreateCommandWidget(aCommand, newParentWidget, isDialog)) {
      return true;
    }
  }
  else {
    CreateVisCommandGroupAndToolBox(aCommand, newParentWidget, aDepthLevel - 1, isDialog);
  }

  return true;
}

/** Create a widget with the command parameters inside
    @param command: command line
    @parent : parent widget
    @isDialog : true if we want apply/cancel button and close at end, false if we want only apply
*/
G4bool G4UIQt::CreateCommandWidget(G4UIcommand* aCommand, QWidget* aParent, G4bool isDialog)
{
  if (aCommand == nullptr) {
    return false;
  }

  // parameters
  auto n_parameterEntry = (G4int)aCommand->GetParameterEntries();
  if (n_parameterEntry > 0) {
    G4UIparameter* param;

    // Re-implementation of G4UIparameter.cc
    auto paramWidget = new QWidget();
    auto gridLayout = new QGridLayout();
    paramWidget->setLayout(gridLayout);

    // Special case for colour, try to display a color chooser if we found red/green/blue parameter
    unsigned int nbColorParameter = 0;
    G4bool isStillColorParameter = false;
    G4bool isColorDialogAdded = false;
    QLabel* redLabel = nullptr;
    QLabel* greenLabel = nullptr;
    QString redDefaultStr = "";
    QString greenDefaultStr = "";
    QString blueDefaultStr = "";
    QWidget* redInput = nullptr;
    QWidget* greenInput = nullptr;

    for (G4int i_thParameter = 0; i_thParameter < n_parameterEntry; i_thParameter++) {
      QString txt;
      param = aCommand->GetParameter(i_thParameter);
      auto label = new QLabel(QString((char*)(param->GetParameterName()).data()));

      if ((label->text() == "red") || (label->text() == "red_or_string")) {
        nbColorParameter++;
        isStillColorParameter = true;
      }
      else if ((label->text() == "green") && isStillColorParameter) {
        nbColorParameter++;
      }
      else if ((label->text() == "blue") && isStillColorParameter) {
        nbColorParameter++;
      }
      else if (! isColorDialogAdded) {
        // not following red/green/blue parameters ?
        if (nbColorParameter == 1) {
          gridLayout->addWidget(redLabel, i_thParameter - 1, 0);
          gridLayout->addWidget(redInput, i_thParameter - 1, 1);
        }
        else if (nbColorParameter == 2) {
          gridLayout->addWidget(redLabel, i_thParameter - 2, 0);
          gridLayout->addWidget(redInput, i_thParameter - 2, 1);
          gridLayout->addWidget(greenLabel, i_thParameter - 1, 0);
          gridLayout->addWidget(greenInput, i_thParameter - 1, 1);
        }
        nbColorParameter = 0;
      }
      // Check parameter type, could be NULL if not found
      QWidget* input = nullptr;
      if ((QString(QChar(param->GetParameterType())) == "d") ||
          (QString(QChar(param->GetParameterType())) == "i"))
      {
        input = new QLineEdit();
        // set default value
        dynamic_cast<QLineEdit*>(input)->setText(QString((char*)(param->GetDefaultValue()).data()));

        if (((label->text() == "red") || (label->text() == "red_or_string")) &&
            isStillColorParameter)
        {
          redDefaultStr = QString((char*)(param->GetDefaultValue()).data());
        }
        else if ((label->text() == "green") && isStillColorParameter) {
          greenDefaultStr = QString((char*)(param->GetDefaultValue()).data());
        }
        else if ((label->text() == "blue") && isStillColorParameter) {
          blueDefaultStr = QString((char*)(param->GetDefaultValue()).data());
        }
      }
      else if (QString(QChar(param->GetParameterType())) == "b") {
        input = new QWidget();
        auto layout = new QHBoxLayout();
        input->setLayout(layout);

        auto buttons = new QButtonGroup();
        auto radioOff = new QRadioButton("0");
        auto radioOn = new QRadioButton("1");
        buttons->addButton(radioOn);
        buttons->addButton(radioOff);
        layout->addWidget(radioOn);
        layout->addWidget(radioOff);

        // set default value
        QString defaultValue = QString((char*)(param->GetDefaultValue()).data());
        if (defaultValue == "0") {
          radioOff->setChecked(true);
        }
        else if (defaultValue == "1") {
          radioOn->setChecked(true);
        }
      }
      else if ((QString(QChar(param->GetParameterType())) == "s") &&
               (! param->GetParameterCandidates().empty()))
      {
        input = new QComboBox();
        QString candidates = QString((char*)(param->GetParameterCandidates()).data());
        QStringList list = candidates.split(" ");

        // add all candidates to widget
        QString defaultValue = QString((char*)(param->GetDefaultValue()).data());
        for (int a = 0; a < list.size(); a++) {
          dynamic_cast<QComboBox*>(input)->addItem(list.at(a));
          if (list.at(a) == defaultValue) {
            dynamic_cast<QComboBox*>(input)->setCurrentIndex(a);
          }
        }
      }
      else if ((QString(QChar(param->GetParameterType())) == "s")) {  // string
        input = new QLineEdit();
        // set default value
        dynamic_cast<QLineEdit*>(input)->setText(QString((char*)(param->GetDefaultValue()).data()));
      }
      else if ((QString(QChar(param->GetParameterType())) == "c")) {  // on/off
        input = new QWidget();
        auto layout = new QHBoxLayout();
        input->setLayout(layout);

        auto buttons = new QButtonGroup();
        auto radioOff = new QRadioButton("off");
        auto radioOn = new QRadioButton("on");
        buttons->addButton(radioOn);
        buttons->addButton(radioOff);
        layout->addWidget(radioOn);
        layout->addWidget(radioOff);

        // set default value
        QString defaultValue = QString((char*)(param->GetDefaultValue()).data());
        if (defaultValue == "off") {
          radioOff->setChecked(true);
        }
        else if (defaultValue == "on") {
          radioOn->setChecked(true);
        }
      }
      else {
        input = new QLineEdit();
        dynamic_cast<QLineEdit*>(input)->setText(QString((char*)(param->GetDefaultValue()).data()));
      }

      txt += "\nParameter : " + QString((char*)(param->GetParameterName()).data()) + "\n";
      if (! param->GetParameterGuidance().empty())
        txt += QString((char*)(param->GetParameterGuidance()).data()) + "\n";

      txt += " Parameter type  : " + QString(QChar(param->GetParameterType())) + "\n";
      if (param->IsOmittable()) {
        txt += " Omittable       : True\n";
      }
      else {
        txt += " Omittable       : False\n";
      }
      if (param->GetCurrentAsDefault()) {
        txt += " Default value   : taken from the current value\n";
      }
      else if (! param->GetDefaultValue().empty()) {
        txt += " Default value   : " + QString((char*)(param->GetDefaultValue()).data()) + "\n";
      }
      if (! param->GetParameterRange().empty()) {
        txt += " Parameter range : " + QString((char*)(param->GetParameterRange()).data()) + "\n";
      }
      if (! param->GetParameterCandidates().empty()) {
        txt +=
          " Candidates      : " + QString((char*)(param->GetParameterCandidates()).data()) + "\n";
      }

      if (isStillColorParameter && (nbColorParameter != 0)) {
        if ((label->text() == "red") || (label->text() == "red_or_string")) {
          redLabel = label;
          redInput = input;
        }
        else if (label->text() == "green") {
          greenLabel = label;
          greenInput = input;
        }
        else if (label->text() == "blue") {
          // we have all, then add a color chooser

          // Create a pixmap with the default color
          QColor qc;
          if ((redDefaultStr != "") && (redDefaultStr != "") && (redDefaultStr != "")) {
            qc.setRgbF(
              redDefaultStr.toDouble(), greenDefaultStr.toDouble(), blueDefaultStr.toDouble());
          }
          QPixmap pixmap = QPixmap(QSize(16, 16));
          pixmap.fill(qc);
          QPainter painter(&pixmap);
          painter.setPen(Qt::black);
          painter.drawRect(0, 0, 15, 15);  // Draw contour

          input = new QPushButton("Change color");
          dynamic_cast<QPushButton*>(input)->setIcon(pixmap);
          dynamic_cast<QPushButton*>(input)->setAccessibleName(
            redDefaultStr + " " + greenDefaultStr + " " + blueDefaultStr);
          label = new QLabel("Choose color");

          // less 1 because we have to add one to the row number
          nbColorParameter--;
          gridLayout->addWidget(label, i_thParameter - nbColorParameter, 0);
          input->setToolTip("Select the current color");
          gridLayout->addWidget(input, i_thParameter - nbColorParameter, 1);

          // Connect pushButton to ColorDialog in callback
          connect(dynamic_cast<QPushButton*>(input), &QPushButton::clicked,
            [this, input]() { this->ChangeColorCallback(input); });
          isColorDialogAdded = true;
          isStillColorParameter = false;
        }
      }
      else {
        gridLayout->addWidget(label, i_thParameter - nbColorParameter, 0);
        input->setToolTip(txt);
        gridLayout->addWidget(input, i_thParameter - nbColorParameter, 1);
      }
    }
    // add command name in hidden value at last line position 0
    auto name = new QLabel(QString((char*)(aCommand->GetCommandPath().data())));
    name->hide();
    gridLayout->addWidget(name, n_parameterEntry - nbColorParameter, 0);

    auto applyButton = new QPushButton("Apply");
    if (! isDialog) {
      gridLayout->addWidget(applyButton, n_parameterEntry - nbColorParameter, 1);
      connect(applyButton, &QPushButton::clicked,
        [this, paramWidget]() { this->VisParameterCallback(paramWidget); });
    }
    else {
      // Apply/Cancel buttons

      applyButton->setAutoDefault(true);
      applyButton->setDefault(true);

      auto cancelButton = new QPushButton(tr("&Cancel"));
      cancelButton->setAutoDefault(true);
      gridLayout->addWidget(cancelButton, n_parameterEntry - nbColorParameter, 1);
      gridLayout->addWidget(applyButton, n_parameterEntry - nbColorParameter, 0);

      connect(applyButton, &QPushButton::clicked,
        [this, paramWidget]() { this->VisParameterCallback(paramWidget); });

      QWidget* parentCheck = aParent;
      QDialog* parentDialog = nullptr;
      G4bool found = false;
      while ((parentCheck->parentWidget()) != nullptr) {
        parentCheck = parentCheck->parentWidget();
        parentDialog = dynamic_cast<QDialog*>(parentCheck);
        if (parentDialog != nullptr) {
          connect(applyButton, SIGNAL(clicked()), parentDialog, SLOT(accept()));
          connect(cancelButton, SIGNAL(clicked()), parentDialog, SLOT(reject()));
          found = true;
        }
      }
      if (! found) {
        return false;
      }
    }

    if (aParent->layout() == nullptr) {
      aParent->setLayout(new QVBoxLayout());
    }
    aParent->layout()->addWidget(paramWidget);
  }

  return true;
}

/** Find a treeItemWidget in the help tree
    @param aCommand item's String to look for
    @return item if found, NULL if not
*/
QTreeWidgetItem* G4UIQt::FindTreeItem(QTreeWidgetItem* aParent, const QString& aCommand)
{
  if (aParent == nullptr) return nullptr;

  // Suppress last "/"
  QString myCommand = aCommand;

  if (myCommand.lastIndexOf("/") == (myCommand.size() - 1)) {
    myCommand = myCommand.left(myCommand.size() - 1);
  }

  if (GetLongCommandPath(aParent) == myCommand) return aParent;

  QTreeWidgetItem* tmp = nullptr;
  for (G4int a = 0; a < aParent->childCount(); ++a) {
    if (tmp == nullptr) tmp = FindTreeItem(aParent->child(a), myCommand);
  }
  return tmp;
}

/**   Build the command list parameters in a QString<br>
 Reimplement partialy the G4UIparameter.cc
 @param aCommand : command to list parameters
 @see G4UIparameter::List()
 @see G4UIcommand::List()
 @return the command list parameters, or "" if nothing
 */
QString G4UIQt::GetCommandList(const G4UIcommand* aCommand)
{
  QString txt = "";
  if (aCommand == nullptr) return txt;

  G4String commandPath = aCommand->GetCommandPath();
  G4String rangeString = aCommand->GetRange();
  auto n_guidanceEntry = (G4int)aCommand->GetGuidanceEntries();
  auto n_parameterEntry = (G4int)aCommand->GetParameterEntries();

  if ((commandPath.empty()) && (rangeString.empty()) && (n_guidanceEntry == 0) &&
      (n_parameterEntry == 0))
  {
    return txt;
  }

  if ((commandPath.length() - 1) != '/') {
    txt += "Command " + QString((char*)(commandPath).data()) + "\n";
  }
  txt += "Guidance :\n";

  for (G4int i_thGuidance = 0; i_thGuidance < n_guidanceEntry; i_thGuidance++) {
    txt += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
  }
  if (! rangeString.empty()) {
    txt += " Range of parameters : " + QString((char*)(rangeString).data()) + "\n";
  }
  if (n_parameterEntry > 0) {
    G4UIparameter* param;

    // Re-implementation of G4UIparameter.cc

    for (G4int i_thParameter = 0; i_thParameter < n_parameterEntry; i_thParameter++) {
      param = aCommand->GetParameter(i_thParameter);
      txt += "\nParameter : " + QString((char*)(param->GetParameterName()).data()) + "\n";
      if (! param->GetParameterGuidance().empty())
        txt += QString((char*)(param->GetParameterGuidance()).data()) + "\n";
      txt += " Parameter type  : " + QString(QChar(param->GetParameterType())) + "\n";
      if (param->IsOmittable()) {
        txt += " Omittable       : True\n";
      }
      else {
        txt += " Omittable       : False\n";
      }
      if (param->GetCurrentAsDefault()) {
        txt += " Default value   : taken from the current value\n";
      }
      else if (! param->GetDefaultValue().empty()) {
        txt += " Default value   : " + QString((char*)(param->GetDefaultValue()).data()) + "\n";
      }
      if (! param->GetParameterRange().empty()) {
        txt += " Parameter range : " + QString((char*)(param->GetParameterRange()).data()) + "\n";
      }
      if (! param->GetParameterCandidates().empty()) {
        txt +=
          " Candidates      : " + QString((char*)(param->GetParameterCandidates()).data()) + "\n";
      }
    }
  }
  return txt;
}

/**   Build the command list parameters in a QString with HTML<br>
   Reimplement partialy the G4UIparameter.cc
   @param aCommand : command to list parameters
   @see G4UIparameter::List()
   @see G4UIcommand::List()
   @return the command list parameters, or "" if nothing
*/
void G4UIQt::updateHelpArea(const G4UIcommand* aCommand)
{
  if (fParameterHelpLabel == nullptr) return;
  if (fParameterHelpTable == nullptr) return;

  fParameterHelpLabel->setTextInteractionFlags(Qt::NoTextInteraction);
  QString txt;
  if (aCommand == nullptr) return;

  G4String commandPath = aCommand->GetCommandPath();
  G4String rangeString = aCommand->GetRange();
  auto n_guidanceEntry = (G4int)aCommand->GetGuidanceEntries();
  auto n_parameterEntry = (G4int)aCommand->GetParameterEntries();

  if ((commandPath.empty()) && (rangeString.empty()) && (n_guidanceEntry == 0) &&
      (n_parameterEntry == 0))
  {
    return;
  }

  if ((commandPath.length() - 1) != '/') {
    txt += "<b>Command </b> " + QString((char*)(commandPath).data()) + "<br />";
  }
  txt += "<b>Guidance :</b> ";
  QString tmpGuidance = "";
  for (G4int i_thGuidance = 0; i_thGuidance < n_guidanceEntry; i_thGuidance++) {
    tmpGuidance = QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data());
    tmpGuidance = tmpGuidance.toHtmlEscaped();
    tmpGuidance.replace("\n", "<br />");
    txt += tmpGuidance + "<br />";
  }
  if (! rangeString.empty()) {
    QString range = QString((char*)(rangeString).data());
    range = range.toHtmlEscaped();
    txt += "<b>Range of parameters : </b> " + range + "<br />";
  }
  else {
    txt += "<br />";
  }
  fParameterHelpLabel->setHtml(txt);

  if (n_parameterEntry > 0) {
    G4UIparameter* param;

    // Re-implementation of G4UIparameter.cc

    fParameterHelpTable->clear();
    fParameterHelpTable->setRowCount(n_parameterEntry);
    fParameterHelpTable->setColumnCount(8);
    fParameterHelpTable->setHorizontalHeaderLabels(
      QStringList() << tr("") << tr("Parameter") << tr("Guidance") << tr("Type") << tr("Ommitable")
                    << tr("Default") << tr("Range") << tr("Candidate"));
    fParameterHelpTable->setColumnWidth(2, 60);

    fParameterHelpTable->verticalHeader()->setVisible(false);
    fParameterHelpTable->setAlternatingRowColors(true);
    fParameterHelpTable->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    fParameterHelpTable->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
    fParameterHelpTable->setWordWrap(true);

    QTableWidgetItem* t = fParameterHelpTable->horizontalHeaderItem(1);
    QFont fnt = t->font();
    G4int size = fnt.pointSize();
    fnt.setPointSize(size - 2);

    for (G4int a = 0; a < n_parameterEntry; a++) {
      param = aCommand->GetParameter(a);
      fParameterHelpTable->setItem(a, 0, new QTableWidgetItem(QString::number(a + 1)));

      fParameterHelpTable->setItem(
        a, 1, new QTableWidgetItem(QString((char*)(param->GetParameterName()).data())));
      if (! param->GetParameterGuidance().empty()) {
        fParameterHelpTable->setItem(
          a, 2, new QTableWidgetItem(QString((char*)(param->GetParameterGuidance()).data())));
      }
      fParameterHelpTable->setItem(
        a, 3, new QTableWidgetItem(QString(QChar(param->GetParameterType()))));

      if (param->IsOmittable()) {
        fParameterHelpTable->setItem(a, 4, new QTableWidgetItem(QString("True")));
      }
      else {
        fParameterHelpTable->setItem(a, 4, new QTableWidgetItem(QString("False")));
      }
      if (param->GetCurrentAsDefault()) {
        fParameterHelpTable->setItem(
          a, 5, new QTableWidgetItem(QString("taken from the current value")));
      }
      else if (! param->GetDefaultValue().empty()) {
        fParameterHelpTable->setItem(
          a, 5, new QTableWidgetItem(QString((char*)(param->GetDefaultValue()).data())));
      }
      if (! param->GetParameterRange().empty()) {
        fParameterHelpTable->setItem(
          a, 6, new QTableWidgetItem(QString((char*)(param->GetParameterRange()).data())));
      }
      if (! param->GetParameterCandidates().empty()) {
        fParameterHelpTable->setItem(
          a, 7, new QTableWidgetItem(QString((char*)(param->GetParameterCandidates()).data())));
      }
      // tooltips
      for (G4int b = 0; b < 8; ++b) {
        QTableWidgetItem* tmp = fParameterHelpTable->item(a, b);
        if (tmp != nullptr) {
          tmp->setToolTip(tmp->text());
          tmp->setFlags(Qt::NoItemFlags);
        }
      }
      fParameterHelpTable->resizeRowToContents(a);
    }
    for (G4int c = 0; c < 8; ++c) {
      if (c != 2) {
        fParameterHelpTable->resizeColumnToContents(c);
      }
    }
    fParameterHelpLabel->setVisible(true);
    fParameterHelpTable->setVisible(true);
  }
}

/**
   Return true if this command takes almost a number (int, double, bool,
   string) as an input
   or a string with a candidate list
 */
G4bool G4UIQt::IsGUICommand(const G4UIcommand* aCommand)
{
  if (aCommand == nullptr) return false;

  auto n_parameterEntry = (G4int)aCommand->GetParameterEntries();

  if (n_parameterEntry > 0) {
    G4UIparameter* param;

    // Re-implementation of G4UIparameter.cc

    for (G4int i_thParameter = 0; i_thParameter < n_parameterEntry; i_thParameter++) {
      param = aCommand->GetParameter(i_thParameter);
      if (QString(QChar(param->GetParameterType())) == "d") {
        return true;
      }
      if (QString(QChar(param->GetParameterType())) == "b") {
        return true;
      }
      if (QString(QChar(param->GetParameterType())) == "i") {
        return true;
      }
      if (QString(QChar(param->GetParameterType())) == "s") {
        return true;
      }
    }
  }
  return false;
}

/**  Implement G4VBasicShell vurtual function
 */
G4bool G4UIQt::GetHelpChoice(G4int&) { return true; }

/**   Event filter method. Every event from QtApplication goes here.<br/>
   We apply a filter only for the Up and Down Arrow press when the QLineEdit<br/>
   is active. If this filter match, Up arrow we give the previous command<br/>
   and Down arrow will give the next if exist.<br/>
   @param obj Emitter of the event
   @param event Kind of event
*/
G4bool G4UIQt::eventFilter(  // Should stay with a minuscule eventFilter because of Qt
  QObject* aObj, QEvent* aEvent)
{
  G4bool tabKeyPress = false;
  G4bool moveCommandCursor = false;
  if (aObj == nullptr) return false;
  if (aEvent == nullptr) return false;

  if (aObj == fHistoryTBTableList) {
    if (aEvent->type() == QEvent::KeyPress) {
      fCommandArea->setFocus();
    }
  }

  if (aObj == fCompleter->popup()) {
    if (aEvent->type() == QEvent::KeyPress) {
      auto e = static_cast<QKeyEvent*>(aEvent);
      if (e->key() == (Qt::Key_Tab)) {
        tabKeyPress = true;
      }
    }
    else if (aEvent->type() == QEvent::Hide) {
      // Store this value
      QString c = fCommandArea->text();
      fLastCompleteCommand = c.left(c.indexOf("<"));
    }
  }

  if (aObj == fCommandArea) {
    if (aEvent->type() == QEvent::KeyPress) {
      auto e = static_cast<QKeyEvent*>(aEvent);
      if ((e->key() == (Qt::Key_Down)) || (e->key() == (Qt::Key_PageDown)) ||
          (e->key() == (Qt::Key_Up)) || (e->key() == (Qt::Key_PageUp)))
      {
        G4int selection = fHistoryTBTableList->currentRow();
        if (fHistoryTBTableList->count() != 0) {
          if (selection == -1) {
            selection = fHistoryTBTableList->count() - 1;
          }
          else {
            if (e->key() == (Qt::Key_Down)) {
              if (selection < (fHistoryTBTableList->count() - 1)) selection++;
            }
            else if (e->key() == (Qt::Key_PageDown)) {
              selection = fHistoryTBTableList->count() - 1;
            }
            else if (e->key() == (Qt::Key_Up)) {
              if (selection > 0) selection--;
            }
            else if (e->key() == (Qt::Key_PageUp)) {
              selection = 0;
            }
          }
          fHistoryTBTableList->clearSelection();
          fHistoryTBTableList->item(selection)->setSelected(true);
          fHistoryTBTableList->setCurrentItem(fHistoryTBTableList->item(selection));
        }
        moveCommandCursor = true;
      }
      else if (e->key() == (Qt::Key_Tab)) {
        tabKeyPress = true;
      }
      else if (((e->modifiers() == Qt::ControlModifier) || (e->modifiers() == Qt::MetaModifier)) &&
               (e->key() == Qt::Key_A))
      {
        fCommandArea->home(false);
        return true;
      }
      else if (((e->modifiers() == Qt::ControlModifier) || (e->modifiers() == Qt::MetaModifier)) &&
               (e->key() == Qt::Key_E))
      {
        fCommandArea->end(false);
        return true;
      }
    }
    else if (aEvent->type() == QEvent::Paint) {
      if (fLastCompleteCommand != "") {
        fCommandArea->setText(fLastCompleteCommand);
        fLastCompleteCommand = "";
      }
    }
  }
  if (tabKeyPress) {
    G4String ss = Complete(fCommandArea->text().toStdString().c_str());
    fCommandArea->setText((char*)(ss.data()));
    fCommandArea->setFocus();
    // do not pass by parent, it will disable widget tab focus !
    return true;
    // L.Garnier : MetaModifier is CTRL for MAC, but I don't want to put a MAC
    // specific #ifdef
  }

  G4bool res = false;
  // change cursor position if needed
  if (moveCommandCursor) {
    fCommandArea->setCursorPosition((int)fCommandArea->text().length());
    fCommandArea->setCursorPosition(4);
  }
  else {
    // pass the event on to the parent class
    res = QObject::eventFilter(aObj, aEvent);
  }
  return res;
}

void G4UIQt::UpdateCommandCompleter()
{
  if (fCommandArea == nullptr) return;

  // remove previous one
  fCommandArea->setCompleter(nullptr);
  if (fCompleter != nullptr) {
    if (fCompleter->popup() != nullptr) {
      fCompleter->popup()->removeEventFilter(this);
    }
  }

  QStandardItemModel* model = CreateCompleterModel("/");
  fCompleter = new QCompleter(model);

  // set all dir visibles in completion
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4UIcommandTree* commandTreeTop = UI->GetTree();
  G4UIcommandTree* aTree = commandTreeTop->FindCommandTree("/");
  if (aTree != nullptr) {
    int Ndir = aTree->GetTreeEntry();
    fCompleter->setMaxVisibleItems(Ndir);
  }
  fCommandArea->setCompleter(fCompleter);
  fCompleter->popup()->installEventFilter(this);
}

QStandardItemModel* G4UIQt::CreateCompleterModel(G4String aCmd)
{
  QList<QStandardItem*> dirModelList;
  QList<QStandardItem*> commandModelList;
  QList<QStandardItem*> subDirModelList;
  QList<QStandardItem*> subCommandModelList;

  G4String strtmp;
  G4int nMatch = 0;

  G4String pName = aCmd;
  G4String remainingPath = aCmd;
  G4String empty = "";
  G4String matchingPath = empty;

  // find the tree
  auto jpre = pName.rfind('/');
  if (jpre != G4String::npos) pName.erase(jpre + 1);
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4UIcommandTree* commandTreeTop = UI->GetTree();
  G4UIcommandTree* aTree = commandTreeTop->FindCommandTree(pName);
  if (aTree != nullptr) {
    G4int Ndir = aTree->GetTreeEntry();
    G4int Ncmd = aTree->GetCommandEntry();

    // directory ...
    for (G4int idir = 1; idir <= Ndir; ++idir) {
      G4String fpdir = aTree->GetTree(idir)->GetPathName();
      // matching test
      if (fpdir.find(remainingPath, 0) == 0) {
        if (nMatch == 0) {
          matchingPath = fpdir;
        }
        else {
          matchingPath = aTree->GetFirstMatchedString(fpdir, matchingPath);
        }
        nMatch++;

        // append to dir model list
        auto item1 = new QStandardItem(fpdir.data());
        QIcon i = QIcon(*fDirIcon);
        item1->setData(1);  // dir
        item1->setIcon(QIcon(*fDirIcon));
        dirModelList.append(item1);

        // Go recursively
        QStandardItemModel* subModel = CreateCompleterModel(fpdir.data());
        for (G4int a = 0; a < subModel->rowCount(); ++a) {
          // copy item (an item could only be part of one model
          auto tempItem = new QStandardItem(subModel->item(a)->text());
          tempItem->setIcon(subModel->item(a)->icon());
          tempItem->setToolTip(subModel->item(a)->toolTip());
          tempItem->setData(subModel->item(a)->data());

          // dir
          if (tempItem->data() == 1) {
            subModel->item(a);
            subDirModelList.append(tempItem);
          }
          // command
          else if (tempItem->data() == 0) {
            subCommandModelList.append(tempItem);
          }
        }
      }
    }

    // command ...
    G4int n_parameterEntry;
    G4String rangeString;
    G4int n_guidanceEntry;
    G4UIcommand* command;
    G4UIparameter* param;
    std::string tooltip;
    G4String params;

    for (G4int icmd = 1; icmd <= Ncmd; ++icmd) {
      tooltip = "";
      params = " ";
      command = aTree->GetCommand(icmd);
      G4String longCommandName = aTree->GetPathName() + command->GetCommandName();
      rangeString = command->GetRange();
      n_guidanceEntry = (G4int)command->GetGuidanceEntries();
      n_parameterEntry = (G4int)command->GetParameterEntries();

      // matching test
      if (longCommandName.find(remainingPath, 0) == 0) {
        if (nMatch == 0) {
          matchingPath = longCommandName + " ";
        }
        else {
          strtmp = longCommandName + " ";
          matchingPath = aTree->GetFirstMatchedString(matchingPath, strtmp);
        }

        // guidance
        for (G4int i_thGuidance = 0; i_thGuidance < n_guidanceEntry; i_thGuidance++) {
          tooltip += std::string((command->GetGuidanceLine(i_thGuidance)).data());
          if (i_thGuidance < n_guidanceEntry - 1) {
            tooltip += "\n";
          }
        }

        // parameters
        for (G4int a = 0; a < n_parameterEntry; a++) {
          param = command->GetParameter(a);
          if (param->IsOmittable()) {
            params += "[<" + param->GetParameterName() + ">] ";
          }
          else {
            params += "<" + param->GetParameterName() + "> ";
          }
        }
        nMatch++;

        // Append to command model list
        auto item = new QStandardItem(G4String(longCommandName + params).data());
        item->setData(0);  // command
        item->setIcon(QIcon(*fCommandIcon));
        item->setToolTip(tooltip.c_str());

        commandModelList.append(item);
      }
    }
  }

  auto model = new QStandardItemModel();
  // initialize the model
  model->setColumnCount(1);

  // concat models
  for (auto a : dirModelList) {
    model->appendRow(a);
  }
  for (auto a : subDirModelList) {
    model->appendRow(a);
  }
  for (auto a : commandModelList) {
    model->appendRow(a);
  }
  for (auto a : subCommandModelList) {
    model->appendRow(a);
  }

  return model;
}

/***************************************************************************/
//
//             SLOTS DEFINITIONS
//
/***************************************************************************/

/**   Called when user give "help" command.
 */
void G4UIQt::ShowHelpCallback() { TerminalHelp(""); }

/**   Called when user click on clear button. Clear the text Output area
 */
void G4UIQt::ClearButtonCallback()
{
  fCoutTBTextArea->clear();
  fG4OutputString.clear();
}

/**   Called when user exit session
 */
void G4UIQt::ExitSession() { SessionTerminate(); }

void G4UIQt::ExitHelp() const {}

/**   Callback call when "click on a menu entry.<br>
   Send the associated command to geant4
*/
void G4UIQt::CommandEnteredCallback()
{
  // split by any new line character
  fCommandArea->setText(fCommandArea->text().trimmed());
#if (QT_VERSION < QT_VERSION_CHECK(5, 15, 0))
  QStringList list = fCommandArea->text().split(QRegularExpression("[\r\n]"),QString::SkipEmptyParts);
#else
  QStringList list = fCommandArea->text().split(QRegularExpression("[\r\n]"),Qt::SkipEmptyParts);
#endif

  // Apply for all commands
  for (G4int a = 0; a < list.size(); ++a) {
    QString txt(list[a].trimmed());
    if (txt != "") {
      fHistoryTBTableList->addItem(txt);
      fHistoryTBTableList->clearSelection();
      fHistoryTBTableList->setCurrentItem(nullptr);
      fCommandArea->setText("");
      G4Qt* interactorManager = G4Qt::getInstance();
      if (interactorManager != nullptr) {
        interactorManager->FlushAndWaitExecution();
      }

      G4String command = txt.toStdString().c_str();
      if (command.substr(0, 4) != "help") {
        ApplyShellCommand(command, exitSession, exitPause);
      }
      else {
        ActivateCommand(command);
      }
    }
  }
  // set the focus to the command line
  fCommandArea->setFocus();

  // Rebuild help tree
  FillHelpTree();

  // Rebuild command completion
  UpdateCommandCompleter();

  if (exitSession) SessionTerminate();
}

/** Callback when the text in the line edit is changed.
 When a newline is inserted, trigger the Activate Command
 on this text end set unchanged the end of the line after the newline.
 */
void G4UIQt::CommandEditedCallback(const QString&)
{
#if (QT_VERSION < QT_VERSION_CHECK(5, 15, 0))
  QStringList list = fCommandArea->text().split(QRegularExpression("[\r\n]"), QString::SkipEmptyParts);
#else
  QStringList list = fCommandArea->text().split(QRegularExpression("[\r\n]"), Qt::SkipEmptyParts);
#endif

  if (list.size() > 1) {  // trigger ActivateCommand
    for (G4int a = 0; a < list.size() - 1; ++a) {
      // set only the first part
      fCommandArea->setText(list[a]);
      // trigger callback
      CommandEnteredCallback();
    }
    // reset unfinished command
    fCommandArea->setText(list[list.size() - 1]);
  }
}

/** Callback when one of the scene/vis parameters has changed
 */
void G4UIQt::VisParameterCallback(QWidget* widget)
{
  if (widget == nullptr) {
    return;
  }

  // Look in all the Grid layout, but only column 1 (0 is the parameter name)
  auto grid = dynamic_cast<QGridLayout*>(widget->layout());
  if (grid == nullptr) {
    return;
  }
  QString command;
  QWidget* name = grid->itemAtPosition(grid->rowCount() - 1, 0)->widget();
  if (dynamic_cast<QLabel*>(name) == nullptr) {
    return;
  }
  command += (dynamic_cast<QLabel*>(name))->text() + " ";

  for (G4int a = 0; a < grid->rowCount() - 1; ++a) {
    QWidget* widgetTmp = grid->itemAtPosition(a, 1)->widget();

    // 4 kind of widgets : QLineEdit / QComboBox / radioButtonsGroup / QPushButton (color chooser)
    if (widgetTmp != nullptr) {
      if (dynamic_cast<QLineEdit*>(widgetTmp) != nullptr) {
        command += (dynamic_cast<QLineEdit*>(widgetTmp))->text() + " ";
      }
      else if (dynamic_cast<QComboBox*>(widgetTmp) != nullptr) {
        command += (dynamic_cast<QComboBox*>(widgetTmp))
                     ->itemText((dynamic_cast<QComboBox*>(widgetTmp))->currentIndex()) +
                   " ";

        // Color chooser
      }
      else if (dynamic_cast<QPushButton*>(widgetTmp) != nullptr) {
        command += widgetTmp->accessibleName() + " ";

        // Check for Button group
      }
      else if (dynamic_cast<QWidget*>(widgetTmp) != nullptr) {
        if (widgetTmp->layout()->count() > 0) {
          if (dynamic_cast<QRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()) != nullptr) {
            QAbstractButton* checked =
              (dynamic_cast<QRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()))
                ->group()
                ->checkedButton();
            if (checked != nullptr) {
              command += (dynamic_cast<QRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()))
                           ->group()
                           ->checkedButton()
                           ->text() +
                         " ";
            }
          }
        }
      }
    }
  }
  if (command != "") {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI != nullptr) {
      UI->ApplyCommand(command.toStdString().c_str());
    }
  }
}

/**   Callback call when "enter" clicked on the command zone.<br>
      If command has no parameters :send the command to geant4
      Else, open a dialog for parameters input
      @param aCommand
*/
void G4UIQt::ButtonCallback(const QString& aCommand)
{
  G4String ss = G4StrUtil::lstrip_copy(G4String(aCommand.toStdString().c_str()));

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return;
  G4UIcommandTree* treeTop = UI->GetTree();

  G4UIcommand* command = treeTop->FindPath(ss);

  if (command != nullptr) {
    // if is GUI, then open a dialog
    if (IsGUICommand(command)) {
      auto menuParameterDialog = new QDialog();

      if (CreateVisCommandGroupAndToolBox(command, menuParameterDialog, 1, true)) {
        menuParameterDialog->setWindowTitle(aCommand);
        menuParameterDialog->setSizePolicy(QSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum));

        // exec this dialog, apply the command automaticaly, and return
        menuParameterDialog->exec();
        return;
      }
      delete menuParameterDialog;
    }
  }

  ApplyShellCommand(ss, exitSession, exitPause);

  // Rebuild help tree
  FillHelpTree();

  if (exitSession) SessionTerminate();
}

/**   This callback is activated when user selected a item in the help tree
 */
void G4UIQt::HelpTreeClicCallback()
{
  QTreeWidgetItem* item = nullptr;
  if (fHelpTreeWidget == nullptr) return;

  QList<QTreeWidgetItem*> list = fHelpTreeWidget->selectedItems();
  if (list.isEmpty()) return;
  item = list.first();
  if (item == nullptr) return;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return;
  G4UIcommandTree* treeTop = UI->GetTree();

  std::string itemText = GetLongCommandPath(item).toStdString();

  // check if it is a command path
  if (item->childCount() > 0) {
    itemText += "/";
  }
  G4UIcommand* command = treeTop->FindPath(itemText.c_str());

  if (command != nullptr) {
    updateHelpArea(command);
  }
  else {  // this is a command
    G4UIcommandTree* path = treeTop->FindCommandTree(itemText.c_str());
    if (path != nullptr) {
      // this is not a command, this is a sub directory
      // We display the Title
      fParameterHelpLabel->setVisible(true);
      fParameterHelpLabel->setText(path->GetTitle().data());
      fParameterHelpTable->setVisible(false);
    }
  }
}

/**   This callback is activated when user double clic on a item in the help tree
 */
void G4UIQt::HelpTreeDoubleClicCallback()
{
  HelpTreeClicCallback();

  QTreeWidgetItem* item = nullptr;
  if (fHelpTreeWidget == nullptr) return;

  QList<QTreeWidgetItem*> list = fHelpTreeWidget->selectedItems();
  if (list.isEmpty()) return;
  item = list.first();
  if (item == nullptr) return;

  fCommandArea->clear();
  fCommandArea->setText(GetLongCommandPath(item));
}

/**   Callback called when user select an old command in the command history<br>
   Give it to the command area.
*/
void G4UIQt::CommandHistoryCallback()
{
  QListWidgetItem* item = nullptr;
  if (fHistoryTBTableList == nullptr) return;

  QList<QListWidgetItem*> list = fHistoryTBTableList->selectedItems();
  if (list.isEmpty()) return;
  item = list.first();
  if (item == nullptr) return;
  fCommandArea->setText(item->text());
}

void G4UIQt::ThreadComboBoxCallback(int) { CoutFilterCallback(""); }

void G4UIQt::CoutFilterCallback(const QString&)
{
  FilterAllOutputTextArea();

  fCoutTBTextArea->repaint();
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(
    fCoutTBTextArea->verticalScrollBar()->maximum());
}

void G4UIQt::SaveOutputCallback()
{
  QString fileName = QFileDialog::getSaveFileName(
    fMainWindow, "Save console output as...", fLastOpenPath, "Save output as...");
  if (fileName != "") {
    QFile data(fileName);
    if (data.open(QFile::WriteOnly | QFile::Truncate)) {
      QTextStream out(&data);
      out << fCoutTBTextArea->toPlainText();
      out.flush();
    }
    data.close();
  }
}

QString G4UIQt::FilterOutput(
  const G4UIOutputString& output, const QString& currentThread, const QString& filter)
{
#ifdef G4MULTITHREADED
  if ((currentThread == "All") || (currentThread == output.fThread)) {
#else
  if (currentThread == "") {
#endif
    if (output.fText.contains(QRegularExpression(filter))) {
      return output.fText;
    }
  }
  return "";
}

void G4UIQt::FilterAllOutputTextArea()
{
  QString currentThread = "";
#ifdef G4MULTITHREADED
  currentThread = fThreadsFilterComboBox->currentText();
  if (currentThread == "Master") {
    currentThread = "";
  }
#endif
  QString filter = fCoutFilter->text();
  G4String previousOutputStream = "";

  QString pref = "";
  QString post = "";

  fCoutTBTextArea->clear();

  for (auto& out : fG4OutputString) {
    if (FilterOutput(out, currentThread, filter) != "") {
      // changing color ?
      if (out.fOutputStream != previousOutputStream) {
        previousOutputStream = out.fOutputStream;
        if (out.fOutputStream == "info") {
          pref = "";
          post = "";
        }
        else if (out.fOutputStream == "warning") {
          pref = "<font color=\"DarkYellow\">";
          post = "</font>";
        }
        else {
          pref = "<font color=\"Red\">";
          post = "</font>";
        }
      }
      fCoutTBTextArea->append(pref + out.fText + post);
    }
  }
}

/**   Callback called when user give a new string to look for<br>
   Display a list of matching commands descriptions. If no string is set,
   will display the complete help tree
*/
void G4UIQt::LookForHelpStringCallback()
{
  fHelpLine->setText(fHelpLine->text().trimmed());
  QString searchText = fHelpLine->text();

  fParameterHelpLabel->setText("");
  fParameterHelpTable->setVisible(false);
  if (searchText == "") {
    // clear old help tree
    fHelpTreeWidget->clear();

    FillHelpTree();

    return;
  }
  OpenHelpTreeOnCommand(searchText);
}

void G4UIQt::OpenHelpTreeOnCommand(const QString& searchText)
{
  // the help tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == nullptr) return;
  G4UIcommandTree* treeTop = UI->GetTree();

  G4int treeSize = treeTop->GetTreeEntry();

  // clear old help tree
  fHelpTreeWidget->clear();

  // look for new items

  int tmp = 0;

#if (QT_VERSION < QT_VERSION_CHECK(5, 15, 0))
  // Before Qt5.15
  QMap<G4int, QString> commandResultMap;
  QMap<G4int, QString> commandChildResultMap;
  for (G4int a = 0; a < treeSize; ++a) {
    G4UIcommand* command = treeTop->FindPath(treeTop->GetTree(a + 1)->GetPathName().data());
    tmp = GetCommandList(command).count(searchText, Qt::CaseInsensitive);
    if (tmp > 0) {
      commandResultMap.insertMulti(
        tmp, QString((char*)(treeTop->GetTree(a + 1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(treeTop->GetTree(a + 1), searchText);
    // insert new childs
    if (! commandChildResultMap.empty()) {
      QMap<int, QString>::const_iterator i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insertMulti(i.key(), i.value());
        i++;
      }
      commandChildResultMap.clear();
    }
  }
#else
  // Qt5.15 and beyond
  QMultiMap<G4int, QString> commandResultMap;
  QMultiMap<G4int, QString> commandChildResultMap;
  for (G4int a = 0; a < treeSize; ++a) {
    G4UIcommand* command = treeTop->FindPath(treeTop->GetTree(a + 1)->GetPathName().data());
    tmp = (int)GetCommandList(command).count(searchText, Qt::CaseInsensitive);
    if (tmp > 0) {
      commandResultMap.insert(tmp, QString((char*)(treeTop->GetTree(a + 1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(treeTop->GetTree(a + 1), searchText);
    // insert new childs
    if (! commandChildResultMap.empty()) {
      auto i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insert(i.key(), i.value());
        ++i;
      }
      commandChildResultMap.clear();
    }
  }
#endif

  // build new help tree
  fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  fHelpTreeWidget->setColumnCount(2);
  QStringList labels;
  labels << QString("Command") << QString("Match");
  fHelpTreeWidget->setHeaderLabels(labels);

  if (commandResultMap.empty()) {
    fParameterHelpLabel->setText("No match found");
    fParameterHelpTable->setVisible(false);
    return;
  }

  auto i = commandResultMap.constEnd();
  i--;
  // 10 maximum progress values
  G4float multValue = 10.0 / (G4float)(i.key());
  QString progressChar = "|";
  QString progressStr = "|";

  QTreeWidgetItem* newItem;
  G4bool end = false;
  while (! end) {
    if (i == commandResultMap.constBegin()) {
      end = true;
    }
    for (G4int a = 0; a < G4int(i.key() * multValue); ++a) {
      progressStr += progressChar;
    }
    newItem = new QTreeWidgetItem();
    QString commandStr = i.value().trimmed();

    if (commandStr.indexOf("/") == 0) {
      commandStr = commandStr.right(commandStr.size() - 1);
    }

    newItem->setText(0, commandStr);
    newItem->setText(1, progressStr);
    fHelpTreeWidget->addTopLevelItem(newItem);
    newItem->setForeground(1, QBrush(Qt::blue));
    progressStr = "|";
    i--;
  }
  fHelpTreeWidget->resizeColumnToContents(0);
  fHelpTreeWidget->sortItems(1, Qt::DescendingOrder);
  //  fHelpTreeWidget->setColumnWidth(1,10);//resizeColumnToContents (1);
}

#if (QT_VERSION < QT_VERSION_CHECK(5, 15, 0))
// Before Qt5.15
QMap<G4int, QString> G4UIQt::LookForHelpStringInChildTree(
  G4UIcommandTree* aCommandTree, const QString& text)
{
  QMap<G4int, QString> commandResultMap;
  if (aCommandTree == NULL) return commandResultMap;
  // Get the Sub directories
  G4int tmp = 0;
  QMap<G4int, QString> commandChildResultMap;
  for (G4int a = 0; a < aCommandTree->GetTreeEntry(); ++a) {
    const G4UIcommand* command = aCommandTree->GetGuidance();
    tmp = GetCommandList(command).count(text, Qt::CaseInsensitive);
    if (tmp > 0) {
      commandResultMap.insertMulti(
        tmp, QString((char*)(aCommandTree->GetTree(a + 1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(aCommandTree->GetTree(a + 1), text);
    if (! commandChildResultMap.empty()) {
      // insert new childs
      QMap<G4int, QString>::const_iterator i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insertMulti(i.key(), i.value());
        ++i;
      }
      commandChildResultMap.clear();
    }
  }
  // Get the Commands
  for (G4int a = 0; a < aCommandTree->GetCommandEntry(); ++a) {
    const G4UIcommand* command = aCommandTree->GetCommand(a + 1);
    tmp = GetCommandList(command).count(text, Qt::CaseInsensitive);
    if (tmp > 0) {
      commandResultMap.insertMulti(
        tmp, QString((char*)(aCommandTree->GetCommand(a + 1)->GetCommandPath()).data()));
    }
  }
  return commandResultMap;
}
#else
// Qt5.15 and beyond
QMultiMap<G4int, QString> G4UIQt::LookForHelpStringInChildTree(
  G4UIcommandTree* aCommandTree, const QString& text)
{
  QMultiMap<G4int, QString> commandResultMap;
  if (aCommandTree == nullptr) return commandResultMap;
  // Get the Sub directories
  G4int tmp = 0;
  QMultiMap<G4int, QString> commandChildResultMap;
  for (G4int a = 0; a < aCommandTree->GetTreeEntry(); ++a) {
    const G4UIcommand* command = aCommandTree->GetGuidance();
    tmp = (int)GetCommandList(command).count(text, Qt::CaseInsensitive);
    if (tmp > 0) {
      commandResultMap.insert(
        tmp, QString((char*)(aCommandTree->GetTree(a + 1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(aCommandTree->GetTree(a + 1), text);
    if (! commandChildResultMap.empty()) {
      // insert new childs
      auto i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insert(i.key(), i.value());
        ++i;
      }
      commandChildResultMap.clear();
    }
  }
  // Get the Commands
  for (G4int a = 0; a < aCommandTree->GetCommandEntry(); ++a) {
    const G4UIcommand* command = aCommandTree->GetCommand(a + 1);
    tmp = (int)GetCommandList(command).count(text, Qt::CaseInsensitive);
    if (tmp > 0) {
      commandResultMap.insert(
        tmp, QString((char*)(aCommandTree->GetCommand(a + 1)->GetCommandPath()).data()));
    }
  }
  return commandResultMap;
}
#endif

QString G4UIQt::GetShortCommandPath(QString commandPath)
{
  if (commandPath.indexOf("/") == 0) {
    commandPath = commandPath.right(commandPath.size() - 1);
  }

  commandPath = commandPath.right(commandPath.size() - commandPath.lastIndexOf("/", -2) - 1);

  if (commandPath.lastIndexOf("/") == (commandPath.size() - 1)) {
    commandPath = commandPath.left(commandPath.size() - 1);
  }

  return commandPath;
}

QString G4UIQt::GetLongCommandPath(QTreeWidgetItem* item)
{
  if (item == nullptr) return "";

  // rebuild path:
  QString itemText = "";
  itemText = item->text(0);

  while (item->parent() != nullptr) {
    itemText = item->parent()->text(0) + "/" + itemText;
    item = item->parent();
  }
  itemText = "/" + itemText;

  return itemText;
}

void G4UIQt::ChangeColorCallback(QWidget* widget)
{
  if (widget == nullptr) {
    return;
  }

  auto button = dynamic_cast<QPushButton*>(widget);
  if (button == nullptr) {
    return;
  }
  QString value = button->accessibleName();

  QColor old;
  old.setRgbF(value.section(" ", 0, 1).toDouble(), value.section(" ", 1, 2).toDouble(),
    value.section(" ", 2, 3).toDouble());
  QColor color =
    QColorDialog::getColor(old, fUITabWidget, "Change color", QColorDialog::ShowAlphaChannel);

  if (color.isValid()) {
    // rebuild the widget icon
    QPixmap pixmap = QPixmap(QSize(16, 16));
    pixmap.fill(color);
    QPainter painter(&pixmap);
    painter.setPen(Qt::black);
    painter.drawRect(0, 0, 15, 15);  // Draw contour

    button->setAccessibleName(QString::number(color.redF()) + " " +
                              QString::number(color.greenF()) + " " +
                              QString::number(color.blueF()) + " ");
    button->setIcon(pixmap);
  }
}

void G4UIQt::ChangeCursorAction(const QString& action)
{
  // Theses actions should be in the app toolbar

  fMoveSelected = true;
  fPickSelected = true;
  fRotateSelected = true;
  fZoomInSelected = true;
  fZoomOutSelected = true;

  if (fToolbarApp == nullptr) return;
  QList<QAction*> list = fToolbarApp->actions();
  for (auto i : list) {
    if (i->data().toString() == action) {
      i->setChecked(true);
      if (i->data().toString() == "pick") {
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/picking true");
        CreatePickInfosDialog();

        fPickInfosDialog->show();
        fPickInfosDialog->raise();
        fPickInfosDialog->activateWindow();
      }
    }
    else if (i->data().toString() == "move") {
      fMoveSelected = false;
      i->setChecked(false);
    }
    else if (i->data().toString() == "pick") {
      fPickSelected = false;
      i->setChecked(false);
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/picking false");
      if (fPickInfosDialog != nullptr) {
        fPickInfosDialog->hide();
      }
    }
    else if (i->data().toString() == "rotate") {
      fRotateSelected = false;
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_in") {
      fZoomInSelected = false;
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_out") {
      fZoomOutSelected = false;
      i->setChecked(false);
    }
  }
  // FIXME : Should connect this to Vis
}

/* A little bit like "void G4OpenGLQtViewer::toggleDrawingAction(int aAction)"
   But for all viewers, not only Qt

FIXME : Should be a feedback when changing viewer !

 */
void G4UIQt::ChangeSurfaceStyle(const QString& action)
{
  // Theses actions should be in the app toolbar

  if (fToolbarApp == nullptr) return;
  QList<QAction*> list = fToolbarApp->actions();
  for (auto i : list) {
    if (i->data().toString() == action) {
      i->setChecked(true);
    }
    else if (i->data().toString() == "hidden_line_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "hidden_line_and_surface_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "solid") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "wireframe") {
      i->setChecked(false);
    }
  }
  // FIXME : Should connect this to Vis

  if (action == "hidden_line_removal") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style w");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 1");
  }
  else if (action == "hidden_line_and_surface_removal") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style s");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 1");
  }
  else if (action == "solid") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style s");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 0");
  }
  else if (action == "wireframe") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style w");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 0");
  }
}

void G4UIQt::OpenIconCallback(const QString& aParam)
{
  QString aCommand = aParam.left(aParam.indexOf(fStringSeparator));
  QString aLabel = aParam.mid(aParam.indexOf(fStringSeparator) + fStringSeparator.length());

  QString nomFich = QFileDialog::getOpenFileName(fMainWindow, aLabel, fLastOpenPath,
    "Macro files (*.mac);;Geant4 files( *.mac *.g4* *.in);;All (*.*)");
  if (nomFich != "") {
    G4UImanager::GetUIpointer()->ApplyCommand(
      (QString(aCommand) + QString(" ") + nomFich).toStdString().c_str());
    QDir dir;
    fLastOpenPath = dir.absoluteFilePath(nomFich);
  }
}

void G4UIQt::SaveIconCallback(const QString& aParam)
{
  QString aCommand = aParam.left(aParam.indexOf(fStringSeparator));
  QString aLabel = aParam.mid(aParam.indexOf(fStringSeparator) + fStringSeparator.length());

  QString nomFich =
    QFileDialog::getSaveFileName(fMainWindow, aLabel, fLastOpenPath, "Macro files (*.mac)");
  if (nomFich != "") {
    G4UImanager::GetUIpointer()->ApplyCommand(
      (QString(aCommand) + QString(" ") + nomFich).toStdString().c_str());
    QDir dir;
    fLastOpenPath = dir.absoluteFilePath(nomFich);
  }
}

void G4UIQt::CreateViewerPropertiesDialog()
{
  if (fViewerPropertiesDialog != nullptr) {
    return;
  }
  fViewerPropertiesDialog = new QDialog();

  fViewerPropertiesDialog->setWindowTitle("Viewer properties");
  fViewerPropertiesDialog->setSizePolicy(
    QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

  if (fViewerPropertiesWidget == nullptr) {
    fViewerPropertiesWidget = new QWidget();
    auto layoutPropertiesWidget = new QVBoxLayout();
    fViewerPropertiesWidget->setLayout(layoutPropertiesWidget);

    CreateEmptyViewerPropertiesWidget();
  }

  auto layoutDialog = new QVBoxLayout();

  layoutDialog->addWidget(fViewerPropertiesWidget);
  layoutDialog->setContentsMargins(0, 0, 0, 0);
  fViewerPropertiesDialog->setLayout(layoutDialog);
}

void G4UIQt::CreatePickInfosDialog()
{
  if (fPickInfosDialog != nullptr) {
    return;
  }
  fPickInfosDialog = new QDialog();

  fPickInfosDialog->setWindowTitle("Pick infos");
  fPickInfosDialog->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

  if (fPickInfosWidget == nullptr) {
    fPickInfosWidget = new QWidget();
    auto layoutPickInfos = new QVBoxLayout();
    fPickInfosWidget->setLayout(layoutPickInfos);

    CreateEmptyPickInfosWidget();
  }

  auto layoutDialog = new QVBoxLayout();

  layoutDialog->addWidget(fPickInfosWidget);
  layoutDialog->setContentsMargins(0, 0, 0, 0);
  fPickInfosDialog->setLayout(layoutDialog);
  fPickInfosDialog->setWindowFlags(Qt::WindowStaysOnTopHint);
}

void G4UIQt::CreateEmptyViewerPropertiesWidget()
{
  if (fViewerPropertiesWidget == nullptr) return;
  if (fViewerPropertiesWidget->layout() == nullptr) return;
  QLayoutItem* wItem;
  if (fViewerPropertiesWidget->layout()->count() != 0) {
    while ((wItem = fViewerPropertiesWidget->layout()->takeAt(0)) != nullptr) {
      delete wItem->widget();
      delete wItem;
    }
  }
  // Add empty one
  auto label = new QLabel("No viewer - Please open a viewer first");
  fViewerPropertiesWidget->layout()->addWidget(label);
  fViewerPropertiesDialog->setWindowTitle("No viewer");
  fViewerPropertiesDialog->setVisible(false);
}

void G4UIQt::CreateEmptyPickInfosWidget()
{
  QLayoutItem* wItem;
  if (fPickInfosWidget->layout()->count() != 0) {
    while ((wItem = fPickInfosWidget->layout()->takeAt(0)) != nullptr) {
      delete wItem->widget();
      delete wItem;
    }
  }
  // Add empty one
  auto label = new QLabel("Click on the object you want to pick");
  fPickInfosWidget->layout()->addWidget(label);
  fPickInfosDialog->setWindowTitle("Nothing to pick");
}

void G4UIQt::ViewerPropertiesIconCallback(int)
{
  CreateViewerPropertiesDialog();

  fViewerPropertiesDialog->show();
  fViewerPropertiesDialog->raise();
  fViewerPropertiesDialog->activateWindow();
}

void G4UIQt::ChangePerspectiveOrtho(const QString& action)
{
  // Theses actions should be in the app toolbar

  if (fToolbarApp == nullptr) return;
  QList<QAction*> list = fToolbarApp->actions();
  QString checked = "";
  for (auto i : list) {
    if (i->data().toString() == action) {
      i->setChecked(true);
      checked = i->data().toString();
    }
    else if (i->data().toString() == "perspective") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "ortho") {
      i->setChecked(false);
    }
  }

  if ((action == "ortho") && (checked == "ortho")) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/projection o");
  }
  else if ((action == "perspective") && (checked == "perspective")) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/projection p");
  }
}

void G4UIQt::SetIconMoveSelected()
{
  // Theses actions should be in the app toolbar
  fMoveSelected = true;
  fRotateSelected = false;
  fPickSelected = false;
  fZoomInSelected = false;
  fZoomOutSelected = false;

  if (fToolbarApp == nullptr) return;
  QList<QAction*> list = fToolbarApp->actions();
  for (auto i : list) {
    if (i->data().toString() == "move") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "rotate") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "pick") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_in") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_out") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconRotateSelected()
{
  // Theses actions should be in the app toolbar
  fRotateSelected = true;
  fMoveSelected = false;
  fPickSelected = false;
  fZoomInSelected = false;
  fZoomOutSelected = false;

  if (fToolbarApp == nullptr) return;
  QList<QAction*> list = fToolbarApp->actions();
  for (auto i : list) {
    if (i->data().toString() == "rotate") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "move") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "pick") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_in") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_out") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconPickSelected()
{
  // Theses actions should be in the app toolbar
  fPickSelected = true;
  fMoveSelected = false;
  fRotateSelected = false;
  fZoomInSelected = false;
  fZoomOutSelected = false;

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "pick") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "move") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "rotate") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_in") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_out") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconZoomInSelected()
{
  // Theses actions should be in the app toolbar
  fZoomInSelected = true;
  fMoveSelected = false;
  fRotateSelected = false;
  fPickSelected = false;
  fZoomOutSelected = false;

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "zoom_in") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "move") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "rotate") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "pick") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_out") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconZoomOutSelected()
{
  // Theses actions should be in the app toolbar
  fZoomOutSelected = true;
  fMoveSelected = false;
  fRotateSelected = false;
  fPickSelected = false;
  fZoomInSelected = false;

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "zoom_out") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "move") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "rotate") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "pick") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "zoom_in") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconSolidSelected()
{
  // Theses actions should be in the app toolbar

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "solid") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "hidden_line_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "hidden_line_and_surface_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "wireframe") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconWireframeSelected()
{
  // Theses actions should be in the app toolbar

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "wireframe") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "hidden_line_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "hidden_line_and_surface_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "solid") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconHLRSelected()
{
  // Theses actions should be in the app toolbar

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "hidden_line_removal") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "solid") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "hidden_line_and_surface_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "wireframe") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconHLHSRSelected()
{
  // Theses actions should be in the app toolbar

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }

  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "hidden_line_and_surface_removal") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "solid") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "hidden_line_removal") {
      i->setChecked(false);
    }
    else if (i->data().toString() == "wireframe") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconPerspectiveSelected()
{
  // Theses actions should be in the app toolbar

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }
  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "perspective") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "ortho") {
      i->setChecked(false);
    }
  }
}

void G4UIQt::SetIconOrthoSelected()
{
  // Theses actions should be in the app toolbar

  QToolBar* bar = fToolbarApp;
  if (! fDefaultIcons) {
    bar = fToolbarUser;
  }

  if (bar == nullptr) return;

  QList<QAction*> list = bar->actions();
  for (auto i : list) {
    if (i->data().toString() == "ortho") {
      i->setChecked(true);
    }
    else if (i->data().toString() == "perspective") {
      i->setChecked(false);
    }
  }
}

G4QTabWidget::G4QTabWidget(QWidget* aParent, G4int sizeX, G4int sizeY)
  : QTabWidget(aParent),
    fTabSelected(false),
    fLastCreated(-1),
    fPreferedSizeX(sizeX + 6)  // margin left+right
    ,
    fPreferedSizeY(sizeY + 58)  // tab label height + margin left+right
{
  setMinimumSize(100, 100);
  QSizePolicy policy = QSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
  setSizePolicy(policy);
}

G4QTabWidget::G4QTabWidget()
  : QTabWidget(), fTabSelected(false), fLastCreated(-1), fPreferedSizeX(0), fPreferedSizeY(0)
{}

G4UIOutputString::G4UIOutputString(QString text, G4String origine, G4String outputStream)
  : fText(text), fThread(origine)
{
  if (! GetOutputList().contains(QString(" ") + outputStream + " ")) {
    fOutputStream = "info";
  }
  else {
    fOutputStream = outputStream;
  }
}


void G4UIQt::TabCloseCallback(int a)
{
  if (fViewerTabWidget == nullptr) return;

  // get the address of the widget
  QWidget* temp = fViewerTabWidget->widget(a);
  // remove the tab
  fViewerTabWidget->removeTab(a);

  // if last QWidget : Add empty string
  G4bool lastTab = true;
  for (G4int c = 0; c < fViewerTabWidget->count(); ++c) {
    if (fViewerTabWidget->tabText(c).contains("viewer")) {
      lastTab = false;
    }
  }

  if (lastTab) {
    CreateEmptyViewerPropertiesWidget();
  }
  // delete the widget
  delete temp;
}

void G4UIQt::ToolBoxActivated(int a)
{
  if (fUITabWidget->widget(a) == fHelpTBWidget) {
    // Rebuild the help tree
    FillHelpTree();
  }
  else if (fUITabWidget->widget(a) == fSceneTreeWidget) {
    fSceneTreeWidget->setVisible(true);
  }
}

void G4QTabWidget::paintEvent(QPaintEvent*)
{
  if (currentWidget() != nullptr) {
    if (isTabSelected()) {
      //      QCoreApplication::sendPostedEvents () ;

      QString text = tabText(currentIndex());

      if (fLastCreated == -1) {
        auto edit = dynamic_cast<QTextEdit*>(currentWidget());
        if (edit == nullptr) {
          QString paramSelect = QString("/vis/viewer/select ") + text;
          G4UImanager* UI = G4UImanager::GetUIpointer();
          if (UI != nullptr) {
            UI->ApplyCommand(paramSelect.toStdString().c_str());
          }
        }
      }
      else {
        fLastCreated = -1;
      }
      setTabSelected(false);
    }
  }
}

G4UIDockWidget::G4UIDockWidget(QString txt) : QDockWidget(txt) {}

void G4UIDockWidget::closeEvent(QCloseEvent* aEvent)
{
  setFloating(false);

  // prevent from closing
  aEvent->ignore();
  // hide them instead
  hide();
}
