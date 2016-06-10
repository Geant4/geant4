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
// $Id: G4UIQt.cc 79182 2014-02-20 09:06:39Z gcosmo $
//
// L. Garnier

#ifdef G4UI_BUILD_QT_SESSION

#include "G4Types.hh"

#include <string.h>

#include "G4UIQt.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommandStatus.hh"

#include "G4Qt.hh"

#include <qapplication.h>
#include <qmessagebox.h>
#include <qlineedit.h>
#include <qwidget.h>
#include <qmenubar.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qsplitter.h>
#include <qscrollbar.h>
#include <qdialog.h>
#include <qevent.h>
#include <qtextedit.h>
#include <qsignalmapper.h>
#include <qtabwidget.h>
#include <qtabbar.h>
#include <qstringlist.h>

#include <qmainwindow.h>
#include <qmenu.h>
#include <qlistwidget.h>
#include <qtreewidget.h>
#include <qgroupbox.h>
#include <qscrollarea.h>
#include <qtoolbox.h>
#include <qradiobutton.h>
#include <qbuttongroup.h>
#include <qcombobox.h>
#include <qsignalmapper.h>
#include <qpainter.h>
#include <qcolordialog.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qdesktopwidget.h>
#include <qmessagebox.h>


#include <stdlib.h>

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
G4UIQt::G4UIQt (
 int argc
,char** argv
)
:fMainWindow(NULL)
,fCommandLabel(NULL)
,fCommandArea(NULL)
,fCoutTBTextArea(NULL)
,fHelpArea(NULL)
,fUITabWidget(NULL)
,fG4cout("")
,fCoutFilter(NULL)
,fHistoryTBTableList(NULL)
,fHelpTreeWidget(NULL)
,fHelpTBWidget(NULL)
,fHistoryTBWidget(NULL)
,fCoutTBWidget(NULL)
,fSceneTreeComponentsTBWidget(NULL)
,fHelpLine(NULL)
,fViewerTabWidget(NULL)
,fCoutText("Output")
,fEmptyViewerTabLabel(NULL)
,fMainSplitterWidget(NULL)
,fRightSplitterWidget(NULL)
,fLeftSplitterWidget(NULL)
,fHelpVSplitter(NULL)
,fViewerTabHandleWidget(NULL)
,fToolbarApp(NULL)
,fToolbarUser(NULL)
,fStringSeparator("__$$$@%%###__")
,fLastOpenPath("")
,fMoveSelected(false)
,fRotateSelected(true)
,fPickSelected(false)
,fZoomInSelected(false)
,fZoomOutSelected(false)
{

  G4Qt* interactorManager = G4Qt::getInstance (argc,argv,(char*)"Qt");
  if (!(QApplication*)interactorManager->GetMainInteractor()) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout        << "G4UIQt : Unable to init Qt. Aborted" << G4endl;
    }
  }
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);
  if(UI!=NULL) UI->SetG4UIWindow(this);

  // Check if already define in external app QMainWindow
  bool found = false;
  Q_FOREACH (QWidget *widget, QApplication::allWidgets()) {
    if ((found== false) && (widget->inherits("QMainWindow"))) {
      found = true;
    }
  }

  if (found) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout        << "G4UIQt : Found an external App with a QMainWindow already defined. Aborted" << G4endl;
    }
    return ;
  }
  fMainWindow = new QMainWindow();

  // the splitter
  fMainSplitterWidget = new QSplitter(Qt::Horizontal);

  fMainSplitterWidget->addWidget(CreateLeftSplitterWidget());
  fMainSplitterWidget->addWidget(CreateRightSplitterWidget());
  
  QSizePolicy policy = QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Preferred);
  policy.setHorizontalStretch(2);
  fRightSplitterWidget->setSizePolicy(policy);
  
  policy = QSizePolicy(QSizePolicy::Expanding,QSizePolicy::Minimum);
  policy.setHorizontalStretch(1);
  fLeftSplitterWidget->setSizePolicy(policy);
  
  fMainWindow->setCentralWidget(fMainSplitterWidget);

  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP

  fMainWindow->setWindowTitle(QFileInfo( QCoreApplication::applicationFilePath() ).fileName()); 
  fMainWindow->move(QPoint(50,50));

  // force the size at be correct at the beggining
  // because the widget is not realized yet, the size of the main window is not up to date. But
  // we need it in order to add some viewer inside
  fMainWindow->resize(fLeftSplitterWidget->width()+fRightSplitterWidget->width()+20,
                      fLeftSplitterWidget->height()+fRightSplitterWidget->height()+20);
  
  // Set not visible until session start
 #if QT_VERSION < 0x040200
  fMainWindow->hide();
 #else
  fMainWindow->setVisible(false);
 #endif
}



G4UIQt::~G4UIQt(
) 
{ 
  G4UImanager* UI = G4UImanager::GetUIpointer();  // TO KEEP
  if(UI!=NULL) {  // TO KEEP
    UI->SetSession(NULL);  // TO KEEP
    UI->SetG4UIWindow(NULL);
    UI->SetCoutDestination(NULL);  // TO KEEP
  }
  
  if (fMainWindow!=NULL) {
    delete fMainWindow;
  }
}

/** Create the History ToolBox Widget
 */
QWidget* G4UIQt::CreateHistoryTBWidget(
) 
{
  fHistoryTBWidget = new QWidget();

  QVBoxLayout *layoutHistoryTB = new QVBoxLayout();
  fHistoryTBTableList = new QListWidget();
  fHistoryTBTableList->setSelectionMode(QAbstractItemView::SingleSelection);
  connect(fHistoryTBTableList, SIGNAL(itemSelectionChanged()), SLOT(CommandHistoryCallback()));
  fHistoryTBTableList->installEventFilter(this);

  layoutHistoryTB->addWidget(fHistoryTBTableList);

  fHistoryTBWidget->setLayout(layoutHistoryTB);
  return fHistoryTBWidget;
}


/** Create the Help ToolBox Widget
 */
QWidget* G4UIQt::CreateHelpTBWidget(
) 
{
  fHelpTBWidget = new QWidget();
  
  QWidget *helpWidget = new QWidget();
  QHBoxLayout *helpLayout = new QHBoxLayout();
  QVBoxLayout *vLayout = new QVBoxLayout();
  fHelpVSplitter = new QSplitter(Qt::Horizontal);
  fHelpLine = new QLineEdit();
  helpLayout->addWidget(new QLabel("Search :"));
  helpLayout->addWidget(fHelpLine);
  connect( fHelpLine, SIGNAL( editingFinished () ), this, SLOT( LookForHelpStringCallback() ) );
  
  // Create Help tree
  FillHelpTree();
  
  fHelpArea = new QTextEdit();
  fHelpArea->setReadOnly(true);
  
  // Set layouts
  
  if (fHelpTreeWidget) {
    fHelpVSplitter->addWidget(fHelpTreeWidget);
  }
  fHelpVSplitter->addWidget(fHelpArea);
  
  vLayout->addWidget(helpWidget);
  vLayout->addWidget(fHelpVSplitter,1);
  
  helpWidget->setLayout(helpLayout);
  fHelpTBWidget->setLayout(vLayout);

  return fHelpTBWidget;
}


/** Create the Cout ToolBox Widget
 */
QWidget* G4UIQt::CreateCoutTBWidget(
) 
{
  fCoutTBWidget = new QGroupBox("Output");

  QVBoxLayout *layoutCoutTB = new QVBoxLayout();

  fCoutTBTextArea = new QTextEdit();

  // set font familly and size
  fCoutTBTextArea->setFontFamily("Courier");
  fCoutTBTextArea->setFontPointSize(12);

  fCoutFilter = new QLineEdit();
  QLabel* coutFilterLabel = new QLabel("Filter : ");
  coutFilterLabel->setToolTip("filter output by...");

  QPushButton *coutTBClearButton = new QPushButton("clear output");
  coutTBClearButton->setToolTip("clear output");
  connect(coutTBClearButton, SIGNAL(clicked()), SLOT(ClearButtonCallback()));
  connect(fCoutFilter, SIGNAL(textEdited ( const QString &)), SLOT(CoutFilterCallback( const QString &)));

  fCoutTBTextArea->setReadOnly(true);

  QWidget* coutButtonWidget = new QWidget();
  QHBoxLayout* layoutCoutTBButtons = new QHBoxLayout();
  layoutCoutTBButtons->addWidget(coutTBClearButton);
  layoutCoutTBButtons->addWidget(coutFilterLabel);
  layoutCoutTBButtons->addWidget(fCoutFilter);
  coutButtonWidget->setLayout(layoutCoutTBButtons);

  // reduce margins
  layoutCoutTBButtons->setContentsMargins(3,3,3,0);

  layoutCoutTB->addWidget(fCoutTBTextArea);
  layoutCoutTB->addWidget(coutButtonWidget);

  fCoutTBWidget->setLayout(layoutCoutTB);

  fCoutTBTextArea->setMinimumSize(100,100);
  
  return fCoutTBWidget;
}


/** Create the VisParameters ToolBox Widget
 */
QWidget* G4UIQt::CreateVisParametersTBWidget(
) 
{
  return NULL;
}


/** Create the VisParameters ToolBox Widget
 */
QWidget* G4UIQt::CreateUITabWidget(
) 
{
  fUITabWidget = new QTabWidget();

  // the right splitter 
  fUITabWidget->addTab(CreateSceneTreeComponentsTBWidget(),"Scene tree");
  fUITabWidget->addTab(CreateHelpTBWidget(),"Help");
  fUITabWidget->addTab(CreateHistoryTBWidget(),"History");
  //  fUITabWidget->setCurrentWidget(fSceneTreeComponentsTBWidget);

  connect(fUITabWidget, SIGNAL(currentChanged(int)), SLOT(ToolBoxActivated(int)));

  return fUITabWidget;
}


QWidget* G4UIQt::CreateSceneTreeComponentsTBWidget(){

  fSceneTreeComponentsTBWidget = new QTabWidget();

#if QT_VERSION < 0x040200
  fSceneTreeComponentsTBWidget->hide();
#else
  fSceneTreeComponentsTBWidget->setVisible(false);
#endif

  return fSceneTreeComponentsTBWidget;
}


QWidget* G4UIQt::CreateLeftSplitterWidget(){

  fLeftSplitterWidget = new QWidget();
  QVBoxLayout * layoutLeftSplitterWidget = new QVBoxLayout();
  layoutLeftSplitterWidget->addWidget(CreateUITabWidget());

  fLeftSplitterWidget->setLayout(layoutLeftSplitterWidget);

  fLeftSplitterWidget->resize(200,200);

  return fLeftSplitterWidget;
}


QWidget* G4UIQt::CreateRightSplitterWidget(){

  fRightSplitterWidget = new QSplitter(Qt::Vertical);

  // Set layouts
  QWidget* commandLineWidget = new QWidget();
  QHBoxLayout *layoutCommandLine = new QHBoxLayout();

  // fill them

  fCommandLabel = new QLabel("");

  fCommandArea = new QLineEdit();
  fCommandArea->installEventFilter(this);
  fCommandArea->activateWindow();

  fCommandArea->setFocusPolicy ( Qt::StrongFocus );
  fCommandArea->setFocus(Qt::TabFocusReason);
  fCommandArea->setToolTip("Apply command");


  layoutCommandLine->addWidget(fCommandLabel);
  layoutCommandLine->addWidget(fCommandArea);

  commandLineWidget->setLayout(layoutCommandLine);

  fEmptyViewerTabLabel = new QLabel("If you want to have a Viewer, please use /vis/open commands.");




  // fill right splitter
    
  //  Create an widget to handle OGL widget and label
  fViewerTabHandleWidget = new QWidget();
  QVBoxLayout * viewerTabHandleLayout = new QVBoxLayout();
  viewerTabHandleLayout->addWidget(fEmptyViewerTabLabel);
  fViewerTabHandleWidget->setLayout(viewerTabHandleLayout);


  fRightSplitterWidget->addWidget(fViewerTabHandleWidget);
  fRightSplitterWidget->addWidget(CreateCoutTBWidget());
  fRightSplitterWidget->addWidget(commandLineWidget);

// set the QGLWidget size policy
  QSizePolicy policy = QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Preferred);
  policy.setVerticalStretch(4);
  fViewerTabHandleWidget->setSizePolicy(policy);
  
  fViewerTabHandleWidget->setMinimumSize(40,40);
  
  commandLineWidget->setMinimumSize(50,50);

  // Connect signal
  connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));
  connect(fCommandArea, SIGNAL(textEdited(const QString &)), SLOT(CommandEditedCallback(const QString &)));

  fRightSplitterWidget->resize(200,200);
  return fRightSplitterWidget;
}


/** Get the ViewerComponents ToolBox Widget
 */
QTabWidget* G4UIQt::GetSceneTreeComponentsTBWidget(
)
{
  return fSceneTreeComponentsTBWidget;
}


/**   Add a new tab widget.
  Create the tab if it was not done
*/
bool G4UIQt::AddTabWidget(
 QWidget* aWidget
,QString name
,int sizeX
,int sizeY
)
{
  // Special case for Qt version between 5.0 and 5.1 on Mac OSX
  // Due to a bug in this Qt version, we can't put a OpenGL Widget inside the QTabWidget.
  // A work around is to put it outside. Returning false will fore the wiewer to put the QGLWidget
  // inside a new QWindow.

#ifdef Q_OS_MAC
  #if QT_VERSION < 0x050100
   #if QT_VERSION >= 0x050000
  QString message = QString(
                            "This Qt version [")+qVersion ()+"] has some issues with the OpenGL viewer.\n"+
  "To prevent problems, you are not allowed to open a Store nor Immediate viewer.\n" +
  "\n" +
  "Please upgrade to Qt version >= 5.1\n";
  
  QMessageBox::warning(fMainWindow, tr("Warning"),
                       tr(message.toStdString().c_str()),
                       QMessageBox::Ok);
  return false;
  #endif
  #endif
#endif
  
  if (fViewerTabWidget == NULL) {
    fViewerTabWidget = new G4QTabWidget(fViewerTabHandleWidget, sizeX, sizeY);
    #if QT_VERSION < 0x040500
#else
    fViewerTabWidget->setTabsClosable (true); 
#endif
    
#if QT_VERSION < 0x040200
#else
    fViewerTabWidget->setUsesScrollButtons (true);
#endif
      
#if QT_VERSION < 0x040500
#else
    connect(fViewerTabWidget,   SIGNAL(tabCloseRequested(int)), this, SLOT(TabCloseCallback(int)));
#endif
    connect(fViewerTabWidget, SIGNAL(currentChanged ( int ) ), SLOT(UpdateTabWidget(int))); 
  }

  if (!aWidget) {
    return false;
  }
// Has to be added before we put it into the fViewerTabHandleWidget widget
  fViewerTabWidget->addTab(aWidget,name);

  // Remove QLabel

  // L.Garnier 26/05/2010 : not exactly the same in qt3. Could cause some
  // troubles
  if (fEmptyViewerTabLabel != NULL) {
    int index = -1;
    index = fViewerTabHandleWidget->layout()->indexOf(fEmptyViewerTabLabel);
    if ( index != -1) {
      fViewerTabHandleWidget->layout()->removeWidget(fEmptyViewerTabLabel);
      delete fEmptyViewerTabLabel;
      fEmptyViewerTabLabel = NULL;
      
      fViewerTabHandleWidget->layout()->addWidget(fViewerTabWidget);
    }
  }



  fViewerTabWidget->setCurrentIndex(fViewerTabWidget->count()-1);

  // Set visible
 #if QT_VERSION < 0x040200
   fViewerTabWidget->setLastTabCreated(fViewerTabWidget->currentIndex());
 #else
   fViewerTabWidget->setLastTabCreated(fViewerTabWidget->currentIndex());
 #endif

  // Problems with resize. The widgets are not realy drawn at this step,
  // then we have to force them on order to check the size.
  // try to computer new size in order not to go out of the screen
  
  // size of current tab
  QSize s = QSize(sizeX,sizeY);
  
  QRect screen = QApplication::desktop()->screenGeometry();
  
  if (fMainWindow->width()-fViewerTabWidget->width()+sizeX > screen.width()) {
    s.setWidth(screen.width()-fMainWindow->width()+fViewerTabWidget->width());
  }
  if (fMainWindow->height()-fViewerTabWidget->height()+sizeY > screen.height()-24) { // 24 is the menuBar height on mac
    s.setHeight(screen.height()-fMainWindow->height()+fViewerTabWidget->height()-24);
  }
  int winWidth = fMainWindow->width();
  int winHeight = fMainWindow->height();
  int oldTabWidth = fViewerTabWidget->width();
  int oldTabHeight = fViewerTabWidget->height();
  int newTabWidth = fViewerTabWidget->sizeHint().width();
  int newTabHeight = fViewerTabWidget->sizeHint().height();

  fViewerTabWidget->setPreferredSize(s);
  fMainWindow->resize(winWidth-oldTabWidth+newTabWidth,
                      winHeight-oldTabHeight+newTabHeight);
  return true;
}


void G4UIQt::UpdateTabWidget(int tabNumber) {
  if ( fViewerTabWidget == NULL) {
    fViewerTabWidget = new G4QTabWidget;
  }
  
  fViewerTabWidget->setCurrentIndex(tabNumber);

  // Send this signal to unblock graphic updates !
  fViewerTabWidget->setTabSelected(false);

 #if QT_VERSION < 0x040200
  fViewerTabWidget->show();
 #else
  fViewerTabWidget->setVisible(true);
 #endif

  // This will send a paintEvent to OGL Viewers
  fViewerTabWidget->setTabSelected(true);
}


/** Send resize event to all tabs
 */
void G4UIQt::ResizeTabWidget( QResizeEvent* e) {
  if ( fViewerTabWidget) {
    for (G4int a=0;a<fViewerTabWidget->count() ;a++) {
      fViewerTabWidget->widget(a)->resize(e->size());
    }
  }
}


/**   Start the Qt main loop
*/
G4UIsession* G4UIQt::SessionStart (
)
{
  G4Qt* interactorManager = G4Qt::getInstance ();
  Prompt("Session :");
  exitSession = false;

  QCoreApplication::sendPostedEvents () ;

  #if QT_VERSION < 0x040200
      fMainWindow->show();
  #else
      fMainWindow->setVisible(true);
  #endif
  interactorManager->DisableSecondaryLoop (); // TO KEEP
  if ((QApplication*)interactorManager->GetMainInteractor())
    ((QApplication*)interactorManager->GetMainInteractor())->exec();

  interactorManager->EnableSecondaryLoop ();
  return this;
}


/**   Display the prompt in the prompt area
   @param aPrompt : string to display as the promt label
   //FIXME : probablement inutile puisque le seul a afficher qq chose d'autre
   que "session" est SecondaryLoop()
*/
void G4UIQt::Prompt (
 G4String aPrompt
)
{
  if (!aPrompt) return;

  fCommandLabel->setText((char*)aPrompt.data());
}



void G4UIQt::SessionTerminate (
)
{
  G4Qt* interactorManager = G4Qt::getInstance ();
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
void G4UIQt::PauseSessionStart (
 const G4String& aState
)
{
  if (!aState) return;

  if(aState=="G4_pause> ") {  // TO KEEP
    SecondaryLoop ("Pause, type continue to exit this state"); // TO KEEP
  } // TO KEEP

  if(aState=="EndOfEvent") { // TO KEEP
    // Picking with feed back in event data Done here !!!
    SecondaryLoop ("End of event, type continue to exit this state"); // TO KEEP
  } // TO KEEP
}



/**
   Begin the secondary loop
   @param a_prompt : label to display as the prompt label
 */
void G4UIQt::SecondaryLoop (
 G4String aPrompt
)
{
  if (!aPrompt) return;

  G4Qt* interactorManager = G4Qt::getInstance (); // TO KEEP ?
  Prompt(aPrompt); // TO KEEP
  exitPause = false; // TO KEEP
  while(1) {
    ((QApplication*)interactorManager)->processEvents(QEventLoop::WaitForMoreEvents);
    if(exitPause==true) break; // TO KEEP
  } // TO KEEP
  Prompt("Session :"); // TO KEEP
}



/**
   Receive a cout from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIQt::ReceiveG4cout (
 const G4String& aString
 )
{
  if (!aString) return 0;
  
  QStringList newStr;
  
  // Add to stringList
  newStr = QStringList(QString((char*)aString.data()).trimmed());
  fG4cout += newStr;

  QStringList result = newStr.filter(fCoutFilter->text());

  if (result.join("").isEmpty()) {
    return 0;
  }
  fCoutTBTextArea->append(result.join(""));
  fCoutTBTextArea->repaint();

  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());

  return 0;
}


/**
   Receive a cerr from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIQt::ReceiveG4cerr (
 const G4String& aString
)
{
  if (!aString) return 0;

  QStringList newStr;

  // Add to stringList
  newStr = QStringList(QString((char*)aString.data()).trimmed());
  fG4cout += newStr;
 
  QStringList result = newStr.filter(fCoutFilter->text());

  // Suppress space, \n,\t,\r...
  if (QString(aString.data()).trimmed() != "") {
    if ((G4StateManager::GetStateManager()->GetCurrentState() == G4State_Abort) ||
        (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Quit )) {
      // In case of Abort or Quit, the useful error message should be in the last error message !
      QMessageBox::critical(fMainWindow, "Error",QString(fLastErrMessage.data())+"\n"+aString.data());
    }
  }
  QColor previousColor = fCoutTBTextArea->textColor();
  fCoutTBTextArea->setTextColor(Qt::red);
  fCoutTBTextArea->append(result.join("\n"));
  fCoutTBTextArea->setTextColor(previousColor);
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
  fCoutTBTextArea->repaint();

  if (QString(aString.data()).trimmed() != "") {
    fLastErrMessage = aString;
  }
  return 0;
}



/**
   Add a new menu to the menu bar
   @param aName name of menu
   @param aLabel label to display
 */
void G4UIQt::AddMenu (
 const char* aName
,const char* aLabel
)
{
  if (aName == NULL) return;
  if (aLabel == NULL) return;

  QMenu *fileMenu = new QMenu(aLabel);
  fMainWindow->menuBar()->addMenu(fileMenu); 

  AddInteractor (aName,(G4Interactor)fileMenu);
}


/**
   Add a new button to a menu
   @param aMenu : parent menu
   @param aLabel : label to display
   @param aCommand : command to execute as a callback
 */
void G4UIQt::AddButton (
 const char* aMenu
,const char* aLabel
,const char* aCommand
)
{
  if(aMenu==NULL) return; // TO KEEP
  if(aLabel==NULL) return; // TO KEEP
  if(aCommand==NULL) return; // TO KEEP

  QMenu *parentTmp = (QMenu*)GetInteractor(aMenu);

  if(parentTmp==NULL) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout << "Menu name " << aMenu<< " does not exist, please define it before using it."<< G4endl;
    }
    return;
  }
  
  // Find the command in the command tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();

  G4String cmd = aCommand;
  G4int cmdEndPos = cmd.find_first_of(" \t");
  if(cmdEndPos!=G4int(std::string::npos)) {
    cmd.erase(cmdEndPos);
  }
  
  if(treeTop->FindPath(cmd) == NULL) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout << "Warning: command '"<< cmd <<"' does not exist, please define it before using it."<< G4endl;
    }
  }

  QSignalMapper *signalMapper = new QSignalMapper(this);
  QAction *action = parentTmp->addAction(aLabel, signalMapper, SLOT(map()));

  connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ButtonCallback(const QString&)));
  signalMapper->setMapping(action, QString(aCommand));
}




/**
 special case for the "open" icon. It will open a file selector and map the return file to the given command.
*/
void G4UIQt::AddIcon(const char* aLabel, const char* aIconFile, const char* aCommand, const char* aFileName){
  if(aLabel==NULL) return; // TO KEEP
  // special case, aCommand could be NULL if aIconFile is not user_icon
  if (aCommand==NULL) {
    if (std::string(aIconFile) == "user_icon") {
      return; // TO KEEP
    }
  }
  QPixmap pix;
  bool userToolBar = false;

  if (std::string(aIconFile) == "user_icon") {
    // try to open a file
    pix = QPixmap(aFileName);
    if (pix.isNull()) {
      G4UImanager* UImanager = G4UImanager::GetUIpointer();
      G4int verbose = UImanager->GetVerboseLevel();
      
      if (verbose >= 2) {
        G4cout << "Warning: file '"<< aFileName <<"' is incorrect or does not exist, this command will not be build"<< G4endl;
      }
      return;
    }
    userToolBar = true; 
  } else if (std::string(aIconFile) == "open") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);

  } else if (std::string(aIconFile) == "save") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "move") {
    const char * const xpm[]={
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
  pix = QPixmap(xpm);

  } else if (std::string(aIconFile) == "rotate") {
    const char * const xpm[]={
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
  pix = QPixmap(xpm);

  } else if (std::string(aIconFile) == "pick") {
    const char * const xpm[]={
        "32 32 2 1",                        
          "       c None",                    
          ".      c #000000",                 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "            .                   ", 
          "            ..                  ", 
          "            ...                 ", 
          "            ....                ", 
          "            .....               ", 
          "            ......              ", 
          "            .......             ", 
          "            .......             ", 
          "            ........            ", 
          "            .....               ", 
          "            ......              ", 
          "            ..  ..              ", 
          "            .   ..              ", 
          "                ...             ", 
          "                 ..             ", 
          "                 ..             ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                ", 
          "                                "}
      ;
  pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "zoom_in") {
    const char * const xpm[]={
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
  pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "zoom_out") {
    const char * const xpm[]={
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
  pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "wireframe") {
    const char * const xpm[]={
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
  pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "solid") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "hidden_line_removal") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "hidden_line_and_surface_removal") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "perspective") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);
  } else if (std::string(aIconFile) == "ortho") {
    const char * const xpm[]={
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
    pix = QPixmap(xpm);
  } else {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout << "Parameter"<< aIconFile <<" not defined"<< G4endl;
    }
    return;
  }
  QToolBar *currentToolbar = NULL;
  if (userToolBar) {
    if (fToolbarUser == NULL) {
      fToolbarUser = new QToolBar();
      fToolbarUser->setIconSize (QSize(20,20));
      fMainWindow->addToolBar(Qt::TopToolBarArea, fToolbarUser);
    }
    currentToolbar = fToolbarUser;
  } else {
    if (fToolbarApp == NULL) {
      fToolbarApp = new QToolBar();
      fToolbarApp->setIconSize (QSize(20,20));
      fMainWindow->addToolBar(Qt::TopToolBarArea, fToolbarApp);
    }
    currentToolbar = fToolbarApp;
  }

  QSignalMapper *signalMapper = new QSignalMapper(this);
  QAction *action = currentToolbar->addAction(pix,aLabel, signalMapper, SLOT(map()));
  

  // special cases :"open"
  if (std::string(aIconFile) == "open") {
    connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(OpenIconCallback(const QString &)));
    QString txt = aCommand + fStringSeparator + aLabel;
    signalMapper->setMapping(action, QString(txt));

  // special cases :"save"
  } else if (std::string(aIconFile) == "save") {
    connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(SaveIconCallback(const QString&)));
    QString txt = aCommand + fStringSeparator + aLabel;
    signalMapper->setMapping(action, QString(txt));

  // special cases : cursor style
  } else if ((std::string(aIconFile) == "move") ||
             (std::string(aIconFile) == "rotate") ||
             (std::string(aIconFile) == "pick") ||
             (std::string(aIconFile) == "zoom_out") ||
             (std::string(aIconFile) == "zoom_in")) {
    action->setCheckable(TRUE);
    action->setChecked(TRUE);
    action->setData(aIconFile);

    connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ChangeCursorStyle(const QString&)));
    signalMapper->setMapping(action, QString(aIconFile));

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
  } else if ((std::string(aIconFile) == "hidden_line_removal") ||
             (std::string(aIconFile) == "hidden_line_and_surface_removal") ||
             (std::string(aIconFile) == "solid") ||
             (std::string(aIconFile) == "wireframe")) {
    action->setCheckable(TRUE);
    action->setChecked(TRUE);
    action->setData(aIconFile);
    connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ChangeSurfaceStyle(const QString&)));
    signalMapper->setMapping(action, QString(aIconFile));

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
  } else if ((std::string(aIconFile) == "perspective") ||
             (std::string(aIconFile) == "ortho")) {
    action->setCheckable(TRUE);
    action->setChecked(TRUE);
    action->setData(aIconFile);
    connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ChangePerspectiveOrtho(const QString&)));
    signalMapper->setMapping(action, QString(aIconFile));

    if (std::string(aIconFile) == "perspective") {
      SetIconPerspectiveSelected();
    }
    if (std::string(aIconFile) == "ortho") {
      SetIconOrthoSelected();
    }

  } else {

    // Find the command in the command tree
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if(UI==NULL) return;
    G4UIcommandTree * treeTop = UI->GetTree();
    if (aCommand != NULL) {
      if(treeTop->FindPath(aCommand) == NULL) {
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        G4int verbose = UImanager->GetVerboseLevel();
        
        if (verbose >= 2) {
          G4cout << "Warning: command '"<< aCommand <<"' does not exist, please define it before using it."<< G4endl;
        }
      }
    }
    
    connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ButtonCallback(const QString&)));
    signalMapper->setMapping(action, QString(aCommand));
  }
}



void G4UIQt::ActivateCommand(
 G4String newCommand
)
{
  if (!fHelpTreeWidget) {
    return;
  }
  // Look for the choosen command "newCommand"
  size_t i = newCommand.index(" ");
  G4String targetCom ="";
  if( i != std::string::npos )
    {
      G4String newValue = newCommand(i+1,newCommand.length()-(i+1));
      newValue.strip(G4String::both);
      targetCom = ModifyToFullPathCommand( newValue );
    }
  if (targetCom != "") {
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

  if (! fHelpTreeWidget ) {
    fHelpTreeWidget = new QTreeWidget();
  }

  // build widget
  fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  QStringList labels;
  labels << QString("Command");
  fHelpTreeWidget->setHeaderLabels(labels);


  connect(fHelpTreeWidget, SIGNAL(itemSelectionChanged ()),this, SLOT(HelpTreeClicCallback()));  
  connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked (QTreeWidgetItem*,int)),this, SLOT(HelpTreeDoubleClicCallback()));  

}
/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIQt::FillHelpTree()
{
  if (! fHelpTreeWidget ) {
    InitHelpTreeAndVisParametersWidget();
  }

  QString searchText = fHelpLine->text();

  if (searchText =="") {
    // clear old help tree
    //    fHelpTreeWidget->clear();
  } else {
    return; 
  }

  if (fHelpArea) {
#if QT_VERSION < 0x040200
    fHelpArea->clear();
#else
    fHelpArea->setText("");
#endif
  }

  if (fHelpLine) {
#if QT_VERSION < 0x040200
    fHelpLine->clear();
#else
    fHelpLine->setText("");
#endif
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();

  G4int treeSize = treeTop->GetTreeEntry();
  QTreeWidgetItem * newItem = NULL;
  QString commandText = "";
  for (int a=0;a<treeSize;a++) {
    // Creating new item
    newItem = NULL;

    commandText = QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).trimmed();

    // if already exist, don't create it !
    for (int b=0;b<fHelpTreeWidget->topLevelItemCount();b++) {
      if (!newItem)
        newItem = FindTreeItem(fHelpTreeWidget->topLevelItem(b),commandText);
    }

    if (newItem == NULL) {
      
      newItem = new QTreeWidgetItem();
      newItem->setText(0,GetShortCommandPath(commandText));
      fHelpTreeWidget->addTopLevelItem(newItem);
    }

    // look for childs
    CreateHelpTree(newItem,treeTop->GetTree(a+1));
  }

}



/**   Fill the Help Tree Widget
   @param aParent : parent item to fill
   @param aCommandTree : commandTree node associate with this part of the Tree
*/
void G4UIQt::CreateHelpTree(
 QTreeWidgetItem *aParent
,G4UIcommandTree *aCommandTree
)
{
  if (aParent == NULL) return;
  if (aCommandTree == NULL) return;


  // Creating new item
  QTreeWidgetItem * newItem;

  QString commandText = "";
  // Get the Sub directories
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {

    commandText = QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).trimmed();
    
    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
      newItem = new QTreeWidgetItem();
      newItem->setText(0,GetShortCommandPath(commandText));
      aParent->addChild(newItem);
    }
    CreateHelpTree(newItem,aCommandTree->GetTree(a+1));
  }

  // Get the Commands

  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    QStringList stringList;
    commandText = QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed();

    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
      newItem = new QTreeWidgetItem();
      newItem->setText(0,GetShortCommandPath(commandText));
      aParent->addChild(newItem);

#if QT_VERSION < 0x040202
      fHelpTreeWidget->setItemExpanded(newItem,false); 
#else
      newItem->setExpanded(false);
#endif
    }
  }
}

 


/**
 Add the following command to the corresponding groupbox
 If depthLevel is 1 : create ToolBox
 If depthLevel is 2 or more : create GroupBox
*/
bool G4UIQt::CreateVisCommandGroupAndToolBox(
 G4UIcommand* aCommand
,QWidget* aParent
,int aDepthLevel
,bool isDialog
)
{
  QString commandText = QString((char*)(aCommand->GetCommandPath().data())).section("/",-aDepthLevel);

  if (commandText == NULL) {
    return false;
  }

  // Look if groupBox is create
  //  QGroupBox* gBoxCommandWidget;
  QWidget* newParentWidget = NULL;
  bool found = false;
  QString commandSection = commandText.left(commandText.indexOf("/"));
  
  if (aDepthLevel == 1) {
    QToolBox* currentParent = dynamic_cast<QToolBox*>(aParent);
    if (currentParent != 0){

      // already exists ?
      for (int a=0; a<currentParent->count(); a++) {
        if (currentParent->itemText(a) == commandSection) {
          found = true;
          newParentWidget = currentParent->widget(a);
        }
      }
    }
    // Not found ? create it
    if (!found) {
      newParentWidget = new QGroupBox();
      newParentWidget->setLayout(new QVBoxLayout());
      if (currentParent != 0){
        currentParent->addItem(newParentWidget,commandSection);
      } else {
        if (!aParent->layout()) {
          aParent->setLayout(new QVBoxLayout());
        }
        aParent->layout()->addWidget(newParentWidget);
      }

      if (commandText.indexOf("/") == -1) {
        
        // Guidance
        QString guidance;
        G4int n_guidanceEntry = aCommand->GetGuidanceEntries();
        for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
          guidance += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
        }
        newParentWidget->setToolTip(guidance);
      }
      
      QScrollArea* sc = dynamic_cast<QScrollArea*>(newParentWidget->parent()->parent());
      if (sc != 0) {
        sc->ensureWidgetVisible(newParentWidget);
        
      }
    }
  } else {

    // try to know if this level is already there
    QGroupBox* currentParent = dynamic_cast<QGroupBox*>(aParent);
    if (currentParent != 0){

      // if depth==2, then we add a [more parameters inside] to the toolBoxItem parent
      // QGroupBox > QWidget > QScrollArea > QToolBox
      if (aDepthLevel == 2){
        QToolBox* parentToolBox = dynamic_cast<QToolBox*>(currentParent->parent()->parent()->parent());
        if (parentToolBox != 0) {
          //          parentToolBox->setItemText(parentToolBox->indexOf(currentParent),"[more parameters inside]");
        }
      }
      for (int a=0; a<aParent->layout()->count(); a++) {
        QGroupBox* gb = dynamic_cast<QGroupBox*>(aParent->layout()->itemAt(a)->widget());
        if (gb != 0) {
          if (gb->title() == commandSection) {
            found = true;
            newParentWidget = gb;
          }
        }
      }
    }
    
    // Not found ? create it
    if (!found) {
      newParentWidget = new QGroupBox();
      newParentWidget->setLayout(new QVBoxLayout());
      if (!aParent->layout()) {
        aParent->setLayout(new QVBoxLayout());
      }
      aParent->layout()->addWidget(newParentWidget);

      // set toolTip
      // Guidance
      QString guidance;
      G4int n_guidanceEntry = aCommand->GetGuidanceEntries();
      for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
        guidance += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
      }
      newParentWidget->setToolTip(guidance);
    }
  }
  
  // fill command groupbox
  if (commandText.indexOf("/") == -1) {
    if (CreateCommandWidget(aCommand, newParentWidget,isDialog)) {
      return true;
    }
  } else {
    CreateVisCommandGroupAndToolBox(aCommand,newParentWidget, aDepthLevel-1,isDialog);
  }

  return true;
}



/** Create a widget with the command parameters inside
    @param command: command line
    @parent : parent widget
    @isDialog : true if we want apply/cancel button and close at end, false if we want only apply
*/
bool G4UIQt::CreateCommandWidget(G4UIcommand* aCommand, QWidget* aParent, bool isDialog) {

  if (aCommand == NULL) {
    return false;
  }


  // parameters
  G4int n_parameterEntry = aCommand->GetParameterEntries();
  if( n_parameterEntry > 0 ) {
    G4UIparameter *param;
      
    // Re-implementation of G4UIparameter.cc
    QWidget* paramWidget = new QWidget();
    QGridLayout* gridLayout = new QGridLayout();
    paramWidget->setLayout(gridLayout);
    
    // Special case for colour, try to display a color chooser if we found red/green/blue parameter
    unsigned int nbColorParameter = 0;
    bool isStillColorParameter = false;
    bool isColorDialogAdded = false;
    QLabel* redLabel = NULL;
    QLabel* greenLabel = NULL;
    QString redDefaultStr = "";
    QString greenDefaultStr = "";
    QString blueDefaultStr = "";
    QWidget* redInput = NULL;
    QWidget* greenInput = NULL;

    for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ ) {
      QString txt;
      param = aCommand->GetParameter(i_thParameter);
      QLabel* label = new QLabel(QString((char*)(param->GetParameterName()).data()));

      if ((label->text() == "red") || (label->text() == "red_or_string")){
        nbColorParameter ++;
        isStillColorParameter = true;
      } else if ((label->text() == "green") && isStillColorParameter) {
        nbColorParameter ++;
      } else if ((label->text() == "blue") && isStillColorParameter) {
        nbColorParameter ++;
      } else if (!isColorDialogAdded) {
          
        // not following red/green/blue parameters ?
        if (nbColorParameter == 1) {
          gridLayout->addWidget(redLabel,i_thParameter-1,0);
          gridLayout->addWidget(redInput,i_thParameter-1,1);
        } else if (nbColorParameter == 2) {
          gridLayout->addWidget(redLabel,i_thParameter-2,0);
          gridLayout->addWidget(redInput,i_thParameter-2,1);
          gridLayout->addWidget(greenLabel,i_thParameter-1,0);
          gridLayout->addWidget(greenInput,i_thParameter-1,1);
        }
        nbColorParameter = 0;
      }
      // Check parameter type, could be NULL if not found
      QWidget* input = NULL;
      if ((QString(QChar(param->GetParameterType())) == "d") || (QString(QChar(param->GetParameterType())) == "i")) {
        input = new QLineEdit();
        // set default value
        dynamic_cast<QLineEdit*>(input)->setText(QString((char*)(param->GetDefaultValue()).data()));

        if (((label->text() == "red") || (label->text() == "red_or_string")) && isStillColorParameter) {
          redDefaultStr = QString((char*)(param->GetDefaultValue()).data());
        } else if ((label->text() == "green") && isStillColorParameter) {
          greenDefaultStr = QString((char*)(param->GetDefaultValue()).data());
        } else if ((label->text() == "green") && isStillColorParameter) {
          blueDefaultStr = QString((char*)(param->GetDefaultValue()).data());
        }

      } else if (QString(QChar(param->GetParameterType())) == "b") {
        input = new QWidget();
        QHBoxLayout* layout = new QHBoxLayout();
        input->setLayout(layout);
        
        QButtonGroup* buttons = new QButtonGroup();
        QRadioButton* radioOff = new QRadioButton("0");
        QRadioButton* radioOn = new QRadioButton("1");
        buttons->addButton(radioOn);
        buttons->addButton(radioOff);
        layout->addWidget(radioOn);
        layout->addWidget(radioOff);

        // set default value
        QString defaultValue = QString((char*)(param->GetDefaultValue()).data());
        if (defaultValue == "0") {
          radioOff->setChecked(true);
        } else if (defaultValue == "1") {
          radioOn->setChecked(true);
        }
      } else if ((QString(QChar(param->GetParameterType())) == "s") && (!param->GetParameterCandidates().isNull())) {
        input = new QComboBox();
        QString candidates = QString((char*)(param->GetParameterCandidates()).data());
        QStringList list = candidates.split (" ");

        // add all candidates to widget
        QString defaultValue = QString((char*)(param->GetDefaultValue()).data());
        for (int a=0; a<list.size(); a++) {
          dynamic_cast<QComboBox*>(input)->addItem(list.at(a));
          if (list.at(a) == defaultValue) {
            dynamic_cast<QComboBox*>(input)->setCurrentIndex(a);
          }
        }

      } else if ((QString(QChar(param->GetParameterType())) == "s")) {  // string
        input = new QLineEdit();
        // set default value
        dynamic_cast<QLineEdit*>(input)->setText(QString((char*)(param->GetDefaultValue()).data()));

      } else if ((QString(QChar(param->GetParameterType())) == "c")) {  // on/off
        input = new QWidget();
        QHBoxLayout* layout = new QHBoxLayout();
        input->setLayout(layout);

        QButtonGroup* buttons = new QButtonGroup();
        QRadioButton* radioOff = new QRadioButton("off");
        QRadioButton* radioOn = new QRadioButton("on");
        buttons->addButton(radioOn);
        buttons->addButton(radioOff);
        layout->addWidget(radioOn);
        layout->addWidget(radioOff);

        // set default value
        QString defaultValue = QString((char*)(param->GetDefaultValue()).data());
        if (defaultValue == "off") {
          radioOff->setChecked(true);
        } else if (defaultValue == "on") {
          radioOn->setChecked(true);
        }

      } else {
        input = new QLineEdit();
        dynamic_cast<QLineEdit*>(input)->setText(QString((char*)(param->GetDefaultValue()).data()));
      }
        
      txt += "\nParameter : " + QString((char*)(param->GetParameterName()).data()) + "\n";
      if( ! param->GetParameterGuidance().isNull() )
        txt += QString((char*)(param->GetParameterGuidance()).data())+ "\n" ;

      txt += " Parameter type  : " + QString(QChar(param->GetParameterType())) + "\n";
      if(param->IsOmittable()){
        txt += " Omittable       : True\n";
      } else {
        txt += " Omittable       : False\n";
      }
      if( param->GetCurrentAsDefault() ) {
        txt += " Default value   : taken from the current value\n";
      } else if( ! param->GetDefaultValue().isNull() ) {
        txt += " Default value   : " + QString((char*)(param->GetDefaultValue()).data())+ "\n";
      }
      if( ! param->GetParameterRange().isNull() ) {
        txt += " Parameter range : " + QString((char*)(param->GetParameterRange()).data())+ "\n";
      }
      if( ! param->GetParameterCandidates().isNull() ) {
        txt += " Candidates      : " + QString((char*)(param->GetParameterCandidates()).data())+ "\n";
      }
        
      if (isStillColorParameter && (nbColorParameter != 0)) {
        if ((label->text() == "red") || (label->text() == "red_or_string")) {
          redLabel = label;
          redInput = input;
        } else if (label->text() == "green") {
          greenLabel = label;
          greenInput = input;
        } else if (label->text() == "blue") {

          // we have all, then add a color chooser

          // Create a pixmap with the default color
          QColor qc;
          if ((redDefaultStr != "") && (redDefaultStr != "") && (redDefaultStr != "")) {
            qc.setRgbF(redDefaultStr.toDouble(),
                       greenDefaultStr.toDouble(),
                       blueDefaultStr.toDouble());
          }
          QPixmap pixmap = QPixmap(QSize(16, 16));
          pixmap.fill (qc);
          QPainter painter(&pixmap);
          painter.setPen(Qt::black);
          painter.drawRect(0,0,15,15); // Draw contour
            
          input = new QPushButton("Change color");
          dynamic_cast<QPushButton*>(input)->setIcon(pixmap);
          dynamic_cast<QPushButton*>(input)->setAccessibleName(redDefaultStr+" "+greenDefaultStr+" "+blueDefaultStr);
          label = new QLabel("Choose color");

          // less 1 because we have to add one to the row number
          nbColorParameter--;
          gridLayout->addWidget(label,i_thParameter-nbColorParameter,0);
          input->setToolTip("Select the current color");
          gridLayout->addWidget(input,i_thParameter-nbColorParameter,1);

          // Connect pushButton to ColorDialog in callback
          QSignalMapper* signalMapper = new QSignalMapper(this);
          signalMapper->setMapping(input,input);
          connect(input, SIGNAL(clicked()), signalMapper, SLOT(map()));
          connect(signalMapper, SIGNAL(mapped(QWidget*)),this, SLOT(ChangeColorCallback(QWidget*)));

          isColorDialogAdded = true;
          isStillColorParameter = false;
        }
      } else {
        gridLayout->addWidget(label,i_thParameter-nbColorParameter,0);
        input->setToolTip(txt);
        gridLayout->addWidget(input,i_thParameter-nbColorParameter,1);
      }
    }
    // add command name in hidden value at last line position 0
    QLabel* name = new QLabel(QString((char*)(aCommand->GetCommandPath().data())));
    name->hide();
    gridLayout->addWidget(name,n_parameterEntry-nbColorParameter,0);

    QPushButton* applyButton = new QPushButton("Apply");
    if (!isDialog) {
      
      gridLayout->addWidget(applyButton,n_parameterEntry-nbColorParameter,1);
      
      QSignalMapper* signalMapper = new QSignalMapper(this);
      signalMapper->setMapping(applyButton, paramWidget);
      connect(applyButton, SIGNAL(clicked()), signalMapper, SLOT(map()));
      connect(signalMapper, SIGNAL(mapped(QWidget*)),this, SLOT(VisParameterCallback(QWidget*)));
    } else {
      // Apply/Cancel buttons
      
      applyButton->setAutoDefault( TRUE );
      applyButton->setDefault( TRUE );
      gridLayout->addWidget(applyButton,n_parameterEntry-nbColorParameter,0);

      QPushButton* cancelButton = new QPushButton( tr( "&Cancel" ));
      cancelButton->setAutoDefault( TRUE );
      gridLayout->addWidget(cancelButton,n_parameterEntry-nbColorParameter,1);
      
      QSignalMapper* signalMapper = new QSignalMapper(this);
      signalMapper->setMapping(applyButton, paramWidget);
      connect(applyButton, SIGNAL(clicked()), signalMapper, SLOT(map()));
      connect(signalMapper, SIGNAL(mapped(QWidget*)),this, SLOT(VisParameterCallback(QWidget*)));

      QWidget * parentCheck = aParent;
      QDialog* parentDialog = NULL;
      bool found = false;
      while ((parentCheck->parentWidget()) != NULL) {
        parentCheck = parentCheck->parentWidget();
        parentDialog = dynamic_cast<QDialog*>(parentCheck);
        if (parentDialog) {
          connect( applyButton, SIGNAL( clicked() ), parentDialog, SLOT( accept() ) );
          connect( cancelButton, SIGNAL( clicked() ), parentDialog, SLOT( reject() ) );
          found = true;
        }
      }
      if (!found) {
        return false;
      }
    }
    
    if (!aParent->layout()) {
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
QTreeWidgetItem* G4UIQt::FindTreeItem(
 QTreeWidgetItem *aParent
,const QString& aCommand
)
{
  if (aParent == NULL) return NULL;

  // Suppress last "/"
  QString myCommand = aCommand;
  
  if (myCommand.lastIndexOf("/") == (myCommand.size()-1)) {
    myCommand = myCommand.left(myCommand.size()-1);
  }

  if (GetLongCommandPath(aParent) == myCommand)
    return aParent;
  
  QTreeWidgetItem * tmp = NULL;
  for (int a=0;a<aParent->childCount();a++) {
    if (!tmp)
      tmp = FindTreeItem(aParent->child(a),myCommand);
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
QString G4UIQt::GetCommandList (
 const G4UIcommand *aCommand
)
{

  QString txt ="";
  if (aCommand == NULL)
    return txt;

  G4String commandPath = aCommand->GetCommandPath();
  G4String rangeString = aCommand->GetRange();
  G4int n_guidanceEntry = aCommand->GetGuidanceEntries();
  G4int n_parameterEntry = aCommand->GetParameterEntries();
  
  if ((commandPath == "") && 
      (rangeString == "") &&
      (n_guidanceEntry == 0) &&
      (n_parameterEntry == 0)) {
    return txt;
  }

  if((commandPath.length()-1)!='/') {
    txt += "Command " + QString((char*)(commandPath).data()) + "\n";
  }
  txt += "Guidance :\n";
  
  for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
    txt += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
  }
  if( ! rangeString.isNull() ) {
    txt += " Range of parameters : " + QString((char*)(rangeString).data()) + "\n";
  }
  if( n_parameterEntry > 0 ) {
    G4UIparameter *param;
    
    // Re-implementation of G4UIparameter.cc
    
    for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ ) {
      param = aCommand->GetParameter(i_thParameter);
      txt += "\nParameter : " + QString((char*)(param->GetParameterName()).data()) + "\n";
      if( ! param->GetParameterGuidance().isNull() )
        txt += QString((char*)(param->GetParameterGuidance()).data())+ "\n" ;
      txt += " Parameter type  : " + QString(QChar(param->GetParameterType())) + "\n";
      if(param->IsOmittable()){
        txt += " Omittable       : True\n";
      } else {
        txt += " Omittable       : False\n";
      }
      if( param->GetCurrentAsDefault() ) {
        txt += " Default value   : taken from the current value\n";
      } else if( ! param->GetDefaultValue().isNull() ) {
        txt += " Default value   : " + QString((char*)(param->GetDefaultValue()).data())+ "\n";
      }
      if( ! param->GetParameterRange().isNull() ) {
        txt += " Parameter range : " + QString((char*)(param->GetParameterRange()).data())+ "\n";
      }
      if( ! param->GetParameterCandidates().isNull() ) {
        txt += " Candidates      : " + QString((char*)(param->GetParameterCandidates()).data())+ "\n";
      }
    }
  }
  return txt;
}


/**
   Return true if this command takes almost a number (int, double, bool, 
   string) as an input
   or a string with a candidate list
 */
G4bool G4UIQt::IsGUICommand(
 const G4UIcommand *aCommand
)
{
  if (aCommand == NULL)
    return false;

  G4int n_parameterEntry = aCommand->GetParameterEntries();
  
  if( n_parameterEntry > 0 ) {
    G4UIparameter *param;
    
    // Re-implementation of G4UIparameter.cc
    
    for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ ) {
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
G4bool G4UIQt::GetHelpChoice(
 G4int&
)
{
  return true;
}


/**   Event filter method. Every event from QtApplication goes here.<br/>
   We apply a filter only for the Up and Down Arrow press when the QLineEdit<br/>
   is active. If this filter match, Up arrow we give the previous command<br/>
   and Down arrow will give the next if exist.<br/>
   @param obj Emitter of the event
   @param event Kind of event
*/
bool G4UIQt::eventFilter( // Should stay with a minuscule eventFilter because of Qt
 QObject *aObj
,QEvent *aEvent
)
{
  bool moveCommandCursor = false;
  if (aObj == NULL) return false;
  if (aEvent == NULL) return false;

  if (aObj == fHistoryTBTableList) {
    if (aEvent->type() == QEvent::KeyPress) {
      fCommandArea->setFocus();
    }
  }
  if (aObj == fCommandArea) {
    if (aEvent->type() == QEvent::KeyPress) {
      QKeyEvent *e = static_cast<QKeyEvent*>(aEvent);
      if ((e->key() == (Qt::Key_Down)) ||
          (e->key() == (Qt::Key_PageDown)) ||
          (e->key() == (Qt::Key_Up)) ||
          (e->key() == (Qt::Key_PageUp))) {
        int selection = fHistoryTBTableList->currentRow();
        if (fHistoryTBTableList->count()) {
          if (selection == -1) {
            selection = fHistoryTBTableList->count()-1;
          } else {
            if (e->key() == (Qt::Key_Down)) {
              if (selection <(fHistoryTBTableList->count()-1))
                selection++;
            } else if (e->key() == (Qt::Key_PageDown)) {
              selection = fHistoryTBTableList->count()-1;
            } else if (e->key() == (Qt::Key_Up)) {
              if (selection >0)
                selection --;
            } else if (e->key() == (Qt::Key_PageUp)) {
              selection = 0;
            }
          }
          fHistoryTBTableList->clearSelection();
#if QT_VERSION < 0x040202
          fHistoryTBTableList->setItemSelected(fHistoryTBTableList->item(selection),true);
#else
          fHistoryTBTableList->item(selection)->setSelected(true);
#endif      
          fHistoryTBTableList->setCurrentItem(fHistoryTBTableList->item(selection));
        }
        moveCommandCursor = true;
      } else if (e->key() == (Qt::Key_Tab)) {
        G4String ss = Complete(fCommandArea->text().toStdString().c_str());
        fCommandArea->setText((char*)(ss.data()));

        // do not pass by parent, it will disable widget tab focus !
        return true;
        // L.Garnier : MetaModifier is CTRL for MAC, but I don't want to put a MAC 
        // specific #ifdef
      } else if (((e->modifiers () == Qt::ControlModifier) || (e->modifiers () == Qt::MetaModifier)) && (e->key() == Qt::Key_A)) {
       fCommandArea->home(false);
       return true;
      } else if (((e->modifiers () == Qt::ControlModifier) || (e->modifiers () == Qt::MetaModifier)) && (e->key() == Qt::Key_E)) {
       fCommandArea->end(false);
       return true;
      }
    }
  }
  bool res= false;
  // change cursor position if needed
  if (moveCommandCursor == true) {
    fCommandArea->setCursorPosition ( fCommandArea->text().length() );
    fCommandArea->setCursorPosition (4);
  } else {
    // pass the event on to the parent class
    res = QObject::eventFilter(aObj, aEvent);
  }
  return res;
}




/***************************************************************************/
//
//             SLOTS DEFINITIONS
//
/***************************************************************************/

/**   Called when user give "help" command.
*/
void G4UIQt::ShowHelpCallback (
)
{
  TerminalHelp("");
}


/**   Called when user click on clear button. Clear the text Output area
*/
void G4UIQt::ClearButtonCallback (
)
{
  fCoutTBTextArea->clear();
  fG4cout.clear();
}

/**   Called when user exit session
*/
void G4UIQt::ExitSession (
)
{
  SessionTerminate();
}

void G4UIQt::ExitHelp(
) const
{
}

/**   Callback call when "click on a menu entry.<br>
   Send the associated command to geant4
*/
void G4UIQt::CommandEnteredCallback (
)
{
  // split by any new line character
  QStringList list = fCommandArea->text().split(QRegExp("[\r\n]"),QString::SkipEmptyParts);

  // Apply for all commands
  for (int a=0; a< list.size(); a++) {
    QString txt (list[a].trimmed());
    if (txt != "") {
      fHistoryTBTableList->addItem(txt);
      fHistoryTBTableList->clearSelection();
      fHistoryTBTableList->setCurrentItem(NULL);
      fCommandArea->setText("");
      G4Qt* interactorManager = G4Qt::getInstance ();
      if (interactorManager) {
        interactorManager->FlushAndWaitExecution();
      }
      
      G4String command = txt.toStdString().c_str();
      if (command(0,4) != "help") {
        ApplyShellCommand (command,exitSession,exitPause);
      } else {
        ActivateCommand(command);
      }
    }
  }
  
  // Rebuild help tree
  FillHelpTree();
  
  if(exitSession==true)
    SessionTerminate();
}


/** Callback when the text in the line edit is changed.
 When a newline is inserted, trigger the Activate Command
 on this text end set unchanged the end of the line after the newline.
 */
void G4UIQt::CommandEditedCallback(const QString &)
{
  QStringList list = fCommandArea->text().split(QRegExp("[\r\n]"),QString::SkipEmptyParts);

  if (list.size() > 1) { // trigger ActivateCommand
    for (int a=0; a<list.size()-1; a++) {
      // set only the first part
      fCommandArea->setText(list[a]);
      // trigger callback
      CommandEnteredCallback();
    }
    // reset unfinished command
    fCommandArea->setText(list[list.size()-1]);
  }
}


/** Callback when one of the scene/vis parameters has changed
 */
void G4UIQt::VisParameterCallback(QWidget* widget){
  if (widget == NULL) {
    return;
  }
  
  // Look in all the Grid layout, but only column 1 (0 is the parameter name)
  QGridLayout* grid = dynamic_cast<QGridLayout*>(widget->layout());
  if (grid == 0) {
    return;
  }
  QString command;
#if QT_VERSION < 0x040400
  QWidget* name = grid->itemAt(grid->columnCount()*(grid->rowCount()-2))->widget();
#else
  QWidget* name = grid->itemAtPosition(grid->rowCount()-1,0)->widget();
#endif
  if (dynamic_cast<QLabel*>(name) == 0) {
    return;
  }
  command += (dynamic_cast<QLabel*>(name))->text()+" ";
  
  for (int a=0;a<grid->rowCount()-1; a++) {
#if QT_VERSION < 0x040400
    QWidget* widgetTmp = grid->itemAt(a*grid->columnCount()+1)->widget();
#else
    QWidget* widgetTmp = grid->itemAtPosition(a,1)->widget();
#endif
    
    // 4 kind of widgets : QLineEdit / QComboBox / radioButtonsGroup / QPushButton (color chooser)
    if (widgetTmp != NULL) {

      if (dynamic_cast<QLineEdit*>(widgetTmp) != 0) {
        command += (dynamic_cast<QLineEdit*>(widgetTmp))->text()+" ";

      } else if (dynamic_cast<QComboBox*>(widgetTmp) != 0){
        command += (dynamic_cast<QComboBox*>(widgetTmp))->itemText((dynamic_cast<QComboBox*>(widgetTmp))->currentIndex())+" ";

        // Color chooser 
      } else if (dynamic_cast<QPushButton*>(widgetTmp) != 0){
        command += widgetTmp->accessibleName()+" ";

        // Check for Button group
      } else if (dynamic_cast<QWidget*>(widgetTmp) != 0){
        if (widgetTmp->layout()->count() > 0){
          if (dynamic_cast<QRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()) != 0) {
            QAbstractButton * checked = (dynamic_cast<QRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()))->group()->checkedButton();
            if (checked != 0) {
              command += (dynamic_cast<QRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()))->group()->checkedButton()->text()+" ";
            }
          }
        }

      }
    }
  }
  if (command != "") {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if(UI != NULL)  {
      UI->ApplyCommand(command.toStdString().c_str());
    }
  }
}


/**   Callback call when "enter" clicked on the command zone.<br>
      If command has no parameters :send the command to geant4
      Else, open a dialog for parameters input
      @param aCommand
*/
void G4UIQt::ButtonCallback (
 const QString& aCommand
)
{
  G4String ss = G4String(aCommand.toStdString().c_str());
  ss = ss.strip(G4String::leading);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();

  G4UIcommand* command = treeTop->FindPath(ss);

  if (command) {
    // if is GUI, then open a dialog
    if (IsGUICommand(command)) {
      QDialog* menuParameterDialog = new QDialog();

      if (CreateVisCommandGroupAndToolBox(command,menuParameterDialog,1,true)) {
        menuParameterDialog->setWindowTitle (aCommand);
        menuParameterDialog->setSizePolicy (QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));

        // exec this dialog, apply the command automaticaly, and return
        menuParameterDialog->exec();
        return;
      }
    }
  }

  ApplyShellCommand(ss,exitSession,exitPause);

  // Rebuild help tree
  FillHelpTree();

  if(exitSession==true) 
    SessionTerminate();
}



/**   This callback is activated when user selected a item in the help tree
*/
void G4UIQt::HelpTreeClicCallback (
)
{
  QTreeWidgetItem* item =  NULL;
  if (!fHelpTreeWidget)
    return ;

  if (!fHelpArea)
    return;
  
  QList<QTreeWidgetItem *> list =fHelpTreeWidget->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
  if (!item)
    return;
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();

  std::string itemText = GetLongCommandPath(item).toStdString();

  // check if it is a command path
  if (item->childCount() > 0) {
    itemText +="/";
  }
  G4UIcommand* command = treeTop->FindPath(itemText.c_str());

  if (command) {
 #if QT_VERSION < 0x040200
    fHelpArea->clear();
    fHelpArea->append(GetCommandList(command));
 #else
    fHelpArea->setText(GetCommandList(command));
 #endif
  } else {  // this is a command
    G4UIcommandTree* path = treeTop->FindCommandTree(itemText.c_str());
    if ( path) {
      // this is not a command, this is a sub directory
      // We display the Title
 #if QT_VERSION < 0x040200
      fHelpArea->clear();
      fHelpArea->append(path->GetTitle().data());
 #else
      fHelpArea->setText(path->GetTitle().data());
 #endif
    }
  }
}
 
/**   This callback is activated when user double clic on a item in the help tree
*/
void G4UIQt::HelpTreeDoubleClicCallback (
)
{
  HelpTreeClicCallback();

  QTreeWidgetItem* item =  NULL;
  if (!fHelpTreeWidget)
    return ;

  if (!fHelpArea)
    return;
  
  QList<QTreeWidgetItem *> list =fHelpTreeWidget->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
  if (!item)
    return;

  fCommandArea->clear();
  fCommandArea->setText(GetLongCommandPath(item));
}


/**   Callback called when user select an old command in the command history<br>
   Give it to the command area.
*/
void G4UIQt::CommandHistoryCallback(
)
{
  QListWidgetItem* item =  NULL;
  if (!fHistoryTBTableList)
    return ;

  
  QList<QListWidgetItem *> list =fHistoryTBTableList->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
  if (!item)
    return;
  fCommandArea->setText(item->text());
}


void G4UIQt::CoutFilterCallback(
const QString & text) {

  QStringList result = fG4cout.filter(text);
  fCoutTBTextArea->setPlainText(result.join("\n"));

  fCoutTBTextArea->repaint();
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());

 }

/**   Callback called when user give a new string to look for<br>
   Display a list of matching commands descriptions. If no string is set,
   will display the complete help tree
*/
void G4UIQt::LookForHelpStringCallback(
)
{
  QString searchText = fHelpLine->text();

#if QT_VERSION < 0x040200
  fHelpArea->clear();
#else
  fHelpArea->setText("");
#endif
  if (searchText =="") {
    // clear old help tree
    fHelpTreeWidget->clear();

    FillHelpTree();

    return;
  } else {
    OpenHelpTreeOnCommand(searchText);
  }
}


void G4UIQt::OpenHelpTreeOnCommand(
 const QString & searchText
)
{
  // the help tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();
  
  G4int treeSize = treeTop->GetTreeEntry();

  // clear old help tree
  fHelpTreeWidget->clear();

  // look for new items

  int tmp = 0;

  QMap<int,QString> commandResultMap;
  QMap<int,QString> commandChildResultMap;

  for (int a=0;a<treeSize;a++) {
    G4UIcommand* command = treeTop->FindPath(treeTop->GetTree(a+1)->GetPathName().data());
    tmp = GetCommandList (command).count(searchText,Qt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(treeTop->GetTree(a+1),searchText);
    // insert new childs
    if (!commandChildResultMap.empty()) {
      QMap<int,QString>::const_iterator i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insertMulti(i.key(),i.value());
        i++;
      }
      commandChildResultMap.clear();
    }
  }

  // build new help tree
  fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  fHelpTreeWidget->setColumnCount(2);
  QStringList labels;
  labels << QString("Command") << QString("Match");
  fHelpTreeWidget->setHeaderLabels(labels);

  if (commandResultMap.empty()) {
#if QT_VERSION < 0x040200
    fHelpArea->clear();
    fHelpArea->append("No match found");
#else
    fHelpArea->setText("No match found");
#endif
    return;
  }

  QMap<int,QString>::const_iterator i = commandResultMap.constEnd();
  i--;
  // 10 maximum progress values
  float multValue = 10.0/(float)(i.key());
  QString progressChar = "|";
  QString progressStr = "|";

  QTreeWidgetItem * newItem;
  bool end = false;
  while (!end) {
    if (i == commandResultMap.constBegin()) {
      end = true;
    }
    for(int a=0;a<int(i.key()*multValue);a++) {
      progressStr += progressChar;
    }
    newItem = new QTreeWidgetItem();
    QString commandStr = i.value().trimmed();

    if (commandStr.indexOf("/") == 0) {
      commandStr = commandStr.right(commandStr.size()-1);
    }
      
    newItem->setText(0,commandStr);
    newItem->setText(1,progressStr);
    fHelpTreeWidget->addTopLevelItem(newItem);
#if QT_VERSION < 0x040200
#else
    newItem->setForeground ( 1, QBrush(Qt::blue) );
#endif
    progressStr = "|";
    i--;
  }
  fHelpTreeWidget->resizeColumnToContents (0);
  fHelpTreeWidget->sortItems(1,Qt::DescendingOrder);
  //  fHelpTreeWidget->setColumnWidth(1,10);//resizeColumnToContents (1);
}




QMap<int,QString> G4UIQt::LookForHelpStringInChildTree(
 G4UIcommandTree *aCommandTree
,const QString & text
 )
{
  QMap<int,QString> commandResultMap;
  if (aCommandTree == NULL) return commandResultMap;
  

  // Get the Sub directories
  int tmp = 0;
  QMap<int,QString> commandChildResultMap;
  
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {
    const G4UIcommand* command = aCommandTree->GetGuidance();
    tmp = GetCommandList (command).count(text,Qt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(aCommandTree->GetTree(a+1),text);
    
    if (!commandChildResultMap.empty()) {
      // insert new childs
      QMap<int,QString>::const_iterator i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insertMulti(i.key(),i.value());
        i++;
      }
      commandChildResultMap.clear();
    }
  }
  // Get the Commands
  
  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    const G4UIcommand* command = aCommandTree->GetCommand(a+1);
    tmp = GetCommandList (command).count(text,Qt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()));
    }
    
  }
  return commandResultMap;
}

  
QString G4UIQt::GetShortCommandPath(
QString commandPath
)
{
  if (commandPath.indexOf("/") == 0) {
    commandPath = commandPath.right(commandPath.size()-1);
  }

  commandPath = commandPath.right(commandPath.size()-commandPath.lastIndexOf("/",-2)-1);
 
 if (commandPath.lastIndexOf("/") == (commandPath.size()-1)) {
    commandPath = commandPath.left(commandPath.size()-1);
 }

 return commandPath;
}


QString G4UIQt::GetLongCommandPath(
 QTreeWidgetItem* item
)
{
  if (item == NULL) return "";

  // rebuild path:
  QString itemText = "";
  itemText = item->text(0);

  while (item->parent() != NULL) {
    itemText = item->parent()->text(0)+"/"+itemText;
    item = item->parent();
  }
  itemText = "/"+itemText;
  
  return itemText;
}


void G4UIQt::ChangeColorCallback(QWidget* widget) {
  if (widget == NULL) {
    return;
  }
  
  QPushButton* button = dynamic_cast<QPushButton*>(widget);
  if (button == 0) {
    return;
  }
  QString value = button->accessibleName();

  QColor old;
  old.setRgbF(value.section(" ",0,1).toDouble(),
              value.section(" ",1,2).toDouble(),
              value.section(" ",2,3).toDouble());
#if QT_VERSION < 0x040500
  bool a;
  QColor color = QColor(QColorDialog::getRgba (old.rgba(),&a,fUITabWidget));
#else
  QColor color = QColorDialog::getColor(old,
                                        fUITabWidget,
                                        "Change color",
					QColorDialog::ShowAlphaChannel);
#endif

  
  if (color.isValid()) {
    // rebuild the widget icon
    QPixmap pixmap = QPixmap(QSize(16, 16));
    pixmap.fill (color);
    QPainter painter(&pixmap);
    painter.setPen(Qt::black);
    painter.drawRect(0,0,15,15); // Draw contour

    button->setAccessibleName(QString::number(color.redF())+" "+
                              QString::number(color.greenF())+" "+
                              QString::number(color.blueF())+" "
                              );
    button->setIcon(pixmap);
    

  }
}


void G4UIQt::ChangeCursorStyle(const QString& action) {

  // Theses actions should be in the app toolbar

  fMoveSelected = true;
  fPickSelected = true;
  fRotateSelected = true;
  fZoomInSelected = true;
  fZoomOutSelected = true;
  
  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == action) {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "move") {
      fMoveSelected = false;
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "pick") {
      fPickSelected = false;
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "rotate") {
      fRotateSelected = false;
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_in") {
      fZoomInSelected = false;
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_out") {
      fZoomOutSelected = false;
      list.at(i)->setChecked(FALSE);
    }
  }
  // FIXME : Should connect this to Vis
}


/* A little bit like "void G4OpenGLQtViewer::toggleDrawingAction(int aAction)"
   But for all viewers, not only Qt

FIXME : Should be a feedback when changing viewer !

 */
void G4UIQt::ChangeSurfaceStyle(const QString& action) {

  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == action) {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "hidden_line_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "hidden_line_and_surface_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "solid") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "wireframe") {
      list.at(i)->setChecked(FALSE);
    }
  }
  // FIXME : Should connect this to Vis

  if (action == "hidden_line_removal") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style w");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 1");

  } else if (action == "hidden_line_and_surface_removal") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style s");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 1");

  } else if (action == "solid") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style s");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 0");

  } else if (action == "wireframe") {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/style w");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/hiddenEdge 0");
  }
}


void G4UIQt::OpenIconCallback(const QString& aParam) {

  QString aCommand = aParam.left(aParam.indexOf(fStringSeparator));
  QString aLabel = aParam.mid(aParam.indexOf(fStringSeparator)+fStringSeparator.length());

  QString nomFich = QFileDialog::getOpenFileName(fMainWindow, aLabel, fLastOpenPath, "Macro files (*.mac)");
  if (nomFich != "") {
    G4UImanager::GetUIpointer()->ApplyCommand((QString(aCommand)+ QString(" ")+ nomFich).toStdString().c_str());
    QDir dir;
    fLastOpenPath = dir.absoluteFilePath(nomFich);
  }
}


void G4UIQt::SaveIconCallback(const QString& aParam) {

  QString aCommand = aParam.left(aParam.indexOf(fStringSeparator));
  QString aLabel = aParam.mid(aParam.indexOf(fStringSeparator)+fStringSeparator.length());

  QString nomFich = QFileDialog::getSaveFileName(fMainWindow, aLabel, fLastOpenPath, "Macro files (*.mac)");
  if (nomFich != "") {
    G4UImanager::GetUIpointer()->ApplyCommand((QString(aCommand)+ QString(" ")+nomFich).toStdString().c_str());
    QDir dir;
    fLastOpenPath = dir.absoluteFilePath(nomFich);
  }
}


  
void G4UIQt::ChangePerspectiveOrtho(const QString& action) {

  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  QString checked = "";
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == action) {
      list.at(i)->setChecked(TRUE);
      checked = list.at(i)->data().toString ();
    } else if (list.at(i)->data().toString () == "persepective") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "ortho") {
      list.at(i)->setChecked(FALSE);
    }
  }

  if ((action == "ortho") && (checked == "ortho")) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/projection o");
  } else if ((action == "perspective") && (checked == "perspective")) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/projection p");
  }
}



void G4UIQt::SetIconMoveSelected() {
  // Theses actions should be in the app toolbar
  fMoveSelected = true;
  fRotateSelected = false;
  fPickSelected = false;
  fZoomInSelected = false;
  fZoomOutSelected = false;

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "move") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "rotate") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "pick") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_in") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_out") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}


void G4UIQt::SetIconRotateSelected() {
  // Theses actions should be in the app toolbar
  fRotateSelected = true;
  fMoveSelected = false;
  fPickSelected = false;
  fZoomInSelected = false;
  fZoomOutSelected = false;

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "rotate") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "move") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "pick") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_in") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_out") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}
 

void G4UIQt::SetIconPickSelected() {
  // Theses actions should be in the app toolbar
  fPickSelected = true;
  fMoveSelected = false;
  fRotateSelected = false;
  fZoomInSelected = false;
  fZoomOutSelected = false;

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "pick") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "move") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "rotate") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_in") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_out") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}
 

void G4UIQt::SetIconZoomInSelected() {
  // Theses actions should be in the app toolbar
  fZoomInSelected = true;
  fMoveSelected = false;
  fRotateSelected = false;
  fPickSelected = false;
  fZoomOutSelected = false;

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "zoom_in") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "move") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "rotate") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "pick") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_out") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}
 

void G4UIQt::SetIconZoomOutSelected() {
  // Theses actions should be in the app toolbar
  fZoomOutSelected = true;
  fMoveSelected = false;
  fRotateSelected = false;
  fPickSelected = false;
  fZoomInSelected = false;

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "zoom_out") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "move") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "rotate") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "pick") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "zoom_in") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}


void G4UIQt::SetIconSolidSelected() {
  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "solid") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "hidden_line_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "hidden_line_and_surface_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "wireframe") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}


void G4UIQt::SetIconWireframeSelected() {
  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "wireframe") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "hidden_line_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "hidden_line_and_surface_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "solid") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}


void G4UIQt::SetIconHLRSelected() {
  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "hidden_line_removal") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "solid") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "hidden_line_and_surface_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "wireframe") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}


void G4UIQt::SetIconHLHSRSelected() {
  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "hidden_line_and_surface_removal") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "solid") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "hidden_line_removal") {
      list.at(i)->setChecked(FALSE);
    } else if (list.at(i)->data().toString () == "wireframe") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}


void G4UIQt::SetIconPerspectiveSelected() {
  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "perspective") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "ortho") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}



void G4UIQt::SetIconOrthoSelected() {
  // Theses actions should be in the app toolbar

  if (fToolbarApp == NULL) return; 
  QList<QAction *> list = fToolbarApp->actions ();
  for (int i = 0; i < list.size(); ++i) {
    if (list.at(i)->data().toString () == "ortho") {
      list.at(i)->setChecked(TRUE);
    } else if (list.at(i)->data().toString () == "perspective") {
      list.at(i)->setChecked(FALSE);
    } 
  }
}



G4QTabWidget::G4QTabWidget(
QWidget* aParent,
int sizeX,
int sizeY
):QTabWidget(aParent)
 ,fTabSelected(false)
 ,fLastCreated(-1)
,fPreferedSizeX(sizeX+6)  // margin left+right
,fPreferedSizeY(sizeY+58)  // tab label height + margin left+right
{
  setMinimumSize(100,100);
  QSizePolicy policy = QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Preferred);
  setSizePolicy(policy);
}

G4QTabWidget::G4QTabWidget(
):QTabWidget()
 ,fTabSelected(false)
 ,fLastCreated(-1)
,fPreferedSizeX(0)
,fPreferedSizeY(0)
{
}


  
#if QT_VERSION < 0x040500
void G4UIQt::TabCloseCallback(int){
#else
void G4UIQt::TabCloseCallback(int a){
#endif
#if QT_VERSION < 0x040500
#else
  if (fViewerTabWidget == NULL) return;
  
  // get the address of the widget
  QWidget* temp = fViewerTabWidget->widget(a);
  // remove the tab
  fViewerTabWidget->removeTab (a);

  // delete the widget
  delete temp;

  // if there is no more tab inside the tab widget
  if (fViewerTabWidget->count() == 0) {
    if (fEmptyViewerTabLabel == NULL) {
      fEmptyViewerTabLabel = new QLabel("If you want to have a Viewer, please use /vis/open commands.");
    }
    
    fViewerTabHandleWidget->layout()->removeWidget(fViewerTabWidget);
    
    fViewerTabHandleWidget->layout()->addWidget(fEmptyViewerTabLabel);
    
    fEmptyViewerTabLabel->show();
  }
#endif
}


void G4UIQt::ToolBoxActivated(int a){
  
  if (fUITabWidget->widget(a) == fHelpTBWidget) {
    // Rebuild the help tree
    FillHelpTree();
  } else if (fUITabWidget->widget(a) == fSceneTreeComponentsTBWidget) {
#if QT_VERSION < 0x040200
    fSceneTreeComponentsTBWidget->show();
#else
    fSceneTreeComponentsTBWidget->setVisible(true);
#endif
  }
}


void G4QTabWidget::paintEvent(
QPaintEvent *
)
{

  if (currentWidget()) {

    if ( isTabSelected()) {

      //      QCoreApplication::sendPostedEvents () ;

      QString text = tabText (currentIndex());

      if (fLastCreated == -1) {
        QString paramSelect = QString("/vis/viewer/select ")+text;
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if(UI != NULL)  {
          UI->ApplyCommand(paramSelect.toStdString().c_str());
        }
      } else {
        fLastCreated = -1;
      }
      setTabSelected(false);
    }
  }
}

#endif
