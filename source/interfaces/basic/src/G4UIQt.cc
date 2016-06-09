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
// $Id: G4UIQt.cc,v 1.53 2010/11/02 15:38:51 lgarnier Exp $
// GEANT4 tag $Name: geant4-09-04 $
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
#include <qlineedit.h>
#include <qwidget.h>
#include <qmenubar.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qtoolbox.h>
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
#if QT_VERSION >= 0x040000
#include <qmenu.h>
#include <qlistwidget.h>
#include <qtreewidget.h>
#else
#include <qaction.h>
#include <qheader.h>
#include <qlistview.h>
#include <qpopupmenu.h>
#include <qwidgetlist.h>
#endif



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
:fHelpArea(NULL)
,fG4cout("")
,fHelpTreeWidget(NULL)
,fHelpTBWidget(NULL)
,fHistoryTBWidget(NULL)
,fCoutTBWidget(NULL)
,fVisParametersTBWidget(NULL)
,fViewComponentsTBWidget(NULL)
,fHelpLine(NULL)
,fTabWidget(NULL)
,fCoutText("Output")
{

  G4Qt* interactorManager = G4Qt::getInstance (argc,argv,(char*)"Qt");
  if (!(QApplication*)interactorManager->GetMainInteractor()) {
    G4cout        << "G4UIQt : Unable to init Qt. Aborted" << G4endl;
  }
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);
  if(UI!=NULL) UI->SetG4UIWindow(this);

  // Check if already define in external app QMainWindow
  bool found = false;
#if QT_VERSION < 0x040000
  // theses lines does nothing exept this one "GLWindow = new QDialog(0..."
  // but if I comment them, it doesn't work...
  QWidgetList  *list = QApplication::allWidgets();
  QWidgetListIt it( *list );         // iterate over the widgets
  QWidget * widget;
  while ( (widget=it.current()) != 0 ) {  // for each widget...
    ++it;
    if ((found== false) && (widget->inherits("QMainWindow"))) {
      found = true;
    }
  }
  delete list;                      // delete the list, not the widgets
#else
  foreach (QWidget *widget, QApplication::allWidgets()) {
    if ((found== false) && (widget->inherits("QMainWindow"))) {
      found = true;
    }
  }
#endif

  if (found) {
    G4cout        << "G4UIQt : Found an external App with a QMainWindow already defined. Aborted" << G4endl;
    return ;
  }
  fMainWindow = new QMainWindow();

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::Initialise after main window creation +++++++++++\n");
#endif

  QWidget *mainWidget = new QWidget(fMainWindow);
#if QT_VERSION < 0x040000
  fMyVSplitter = new QSplitter(Qt::Horizontal,fMainWindow);
  fToolBox = new QToolBox(fMyVSplitter);
#else
  fMyVSplitter = new QSplitter(Qt::Horizontal,fMainWindow);
  fToolBox = new QToolBox();
#endif

  // Set layouts

  QWidget* commandLineWidget = new QWidget(mainWidget);
#if QT_VERSION < 0x040000
  QVBoxLayout *layoutCommandLine = new QVBoxLayout(commandLineWidget);
#else
  QVBoxLayout *layoutCommandLine = new QVBoxLayout();
#endif

  // fill them

  fCommandLabel = new QLabel("",commandLineWidget);

  fCommandArea = new QLineEdit(commandLineWidget);
  fCommandArea->installEventFilter(this);
#if QT_VERSION < 0x040000
  fCommandArea->setActiveWindow();
#else
  fCommandArea->activateWindow();
#endif

#if QT_VERSION < 0x040000
  fCommandArea->setFocusPolicy ( QWidget::StrongFocus );
  fCommandArea->setFocus();
#else
  fCommandArea->setFocusPolicy ( Qt::StrongFocus );
  fCommandArea->setFocus(Qt::TabFocusReason);
#endif



  layoutCommandLine->addWidget(fCommandLabel);
  layoutCommandLine->addWidget(fCommandArea);
  QVBoxLayout *mainLayout;
#if QT_VERSION >= 0x040000
  mainLayout = new QVBoxLayout();
#else
  mainLayout = new QVBoxLayout(mainWidget);
#endif

  fHelpTBWidget = new QWidget(fToolBox);
  fHistoryTBWidget = new QWidget(fToolBox);
  fCoutTBWidget = new QWidget(fToolBox);
  fVisParametersTBWidget = new QWidget(fToolBox);
  fViewComponentsTBWidget = new QWidget(fToolBox);
  
  CreateVisParametersTBWidget();
  CreateViewComponentsTBWidget();
  CreateHelpTBWidget();
  CreateCoutTBWidget();
  CreateHistoryTBWidget();

  // the splitter 
  //  fToolBox->addItem(fVisParametersTBWidget,"Vis parameters");
  //  fToolBox->addItem(fViewComponentsTBWidget,"Viewer components");
  fToolBox->addItem(fHelpTBWidget,"Help");
  fToolBox->addItem(fCoutTBWidget,"Cout");
  fToolBox->addItem(fHistoryTBWidget,"History");



  fToolBox->setSizePolicy (QSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed));

#if QT_VERSION < 0x040000
  fEmptyViewerTabLabel = new QLabel(fToolBox,"         If you want to have a Viewer, please use /vis/open commands. ");
#else
  fEmptyViewerTabLabel = new QLabel("         If you want to have a Viewer, please use /vis/open commands. ");
#endif

  // Only at creation. Will be set visible when sessionStart();
#if QT_VERSION >= 0x040000
 #if QT_VERSION >= 0x040200
  fEmptyViewerTabLabel->setVisible(false);
 #else
  fEmptyViewerTabLabel->hide();
 #endif
#else
  fEmptyViewerTabLabel->hide();
#endif


#if QT_VERSION >= 0x040000
  fMyVSplitter->addWidget(fToolBox);
  fMyVSplitter->addWidget(fEmptyViewerTabLabel);
#endif


#if QT_VERSION >= 0x040000
  commandLineWidget->setLayout(layoutCommandLine);
#endif
  commandLineWidget->setSizePolicy (QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));
  mainLayout->addWidget(fMyVSplitter,1);
  mainLayout->addWidget(commandLineWidget);

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt :: 5\n");
#endif

#if QT_VERSION >= 0x040000
  mainWidget->setLayout(mainLayout);
#endif

  fMainWindow->setCentralWidget(mainWidget);

#if QT_VERSION < 0x040000

  // Add a quit subMenu
  QPopupMenu *fileMenu = new QPopupMenu( fMainWindow);
  fileMenu->insertItem( "&Quit",  this, SLOT(ExitSession()), CTRL+Key_Q );
  fMainWindow->menuBar()->insertItem( QString("&File"), fileMenu );

#else

  // Add a quit subMenu
  QMenu *fileMenu = fMainWindow->menuBar()->addMenu("File");
  fileMenu->addAction("Quit", this, SLOT(ExitSession()));

#endif

  AddInteractor ("file",(G4Interactor)fileMenu);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt :: 6\n");
#endif

  // Connect signal
  connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));
  connect(fToolBox, SIGNAL(currentChanged(int)), SLOT(ToolBoxActivated(int)));

  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP

#if QT_VERSION < 0x040000
  fMainWindow->setCaption( tr( "G4UI Session" ));
  fMainWindow->resize(900,600); 
  fMainWindow->move(50,100);
#else
  fMainWindow->setWindowTitle( tr("G4UI Session") ); 
  fMainWindow->resize(900,600); 
  fMainWindow->move(QPoint(50,100));
#endif

  // Set not visible until session start
#if QT_VERSION >= 0x040000
 #if QT_VERSION >= 0x040200
  fMainWindow->setVisible(false);
 #else
  fMainWindow->hide();
 #endif
#else
  fMainWindow->hide();
#endif

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt END\n");
#endif
}



G4UIQt::~G4UIQt(
) 
{ 
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::~G4UIQt Delete\n");
#endif
  G4UImanager* UI = G4UImanager::GetUIpointer();  // TO KEEP
  if(UI!=NULL) {  // TO KEEP
    UI->SetSession(NULL);  // TO KEEP
    UI->SetG4UIWindow(NULL);
    UI->SetCoutDestination(NULL);  // TO KEEP
  }
  
  if (fMainWindow!=NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::~G4UIQt DELETE fMainWindow\n");
#endif
    delete fMainWindow;
  }
}

/** Create the History ToolBox Widget
 */
void G4UIQt::CreateHistoryTBWidget(
) 
{

#if QT_VERSION < 0x040000
  QVBoxLayout *layoutHistoryTB = new QVBoxLayout(fHistoryTBWidget);

  fHistoryTBTableList = new QListView(fHistoryTBWidget);
  fHistoryTBTableList->setSorting (-1, FALSE);
  fHistoryTBTableList->setSelectionMode(QListView::Single);
  fHistoryTBTableList->addColumn("");
  fHistoryTBTableList->header()->hide();
  connect(fHistoryTBTableList, SIGNAL(selectionChanged()), SLOT(CommandHistoryCallback()));
#else
  QVBoxLayout *layoutHistoryTB = new QVBoxLayout();
  fHistoryTBTableList = new QListWidget();
  fHistoryTBTableList->setSelectionMode(QAbstractItemView::SingleSelection);
  connect(fHistoryTBTableList, SIGNAL(itemSelectionChanged()), SLOT(CommandHistoryCallback()));
#endif
  fHistoryTBTableList->installEventFilter(this);

  layoutHistoryTB->addWidget(fHistoryTBTableList);

#if QT_VERSION >= 0x040000
  fHistoryTBWidget->setLayout(layoutHistoryTB);
#endif
}

/** Create the Help ToolBox Widget
 */
void G4UIQt::CreateHelpTBWidget(
) 
{

  
#if QT_VERSION < 0x040000
  QWidget *helpWidget = new QWidget(fHelpTBWidget);
  QHBoxLayout *helpLayout = new QHBoxLayout(helpWidget);
  fHelpVSplitter = new QSplitter(Qt::Horizontal,fHelpTBWidget);
#else
  QWidget *helpWidget = new QWidget();
  QHBoxLayout *helpLayout = new QHBoxLayout();
  QVBoxLayout *vLayout = new QVBoxLayout();
  fHelpVSplitter = new QSplitter(Qt::Horizontal);
#endif
  fHelpLine = new QLineEdit(fHelpTBWidget);
  helpLayout->addWidget(new QLabel("Search :",helpWidget));
  helpLayout->addWidget(fHelpLine);
#if QT_VERSION < 0x040000
  connect( fHelpLine, SIGNAL( returnPressed () ), this, SLOT( LookForHelpStringCallback() ) );
#else
  connect( fHelpLine, SIGNAL( editingFinished () ), this, SLOT( LookForHelpStringCallback() ) );
#endif
  
  // Create Help tree
  FillHelpTree();
  
  fHelpArea = new QTextEdit(fHelpVSplitter);
  fHelpArea->setReadOnly(true);
  
  // Set layouts
  
#if QT_VERSION >= 0x040000
  if (fHelpTreeWidget) {
    fHelpVSplitter->addWidget(fHelpTreeWidget);
  }
  fHelpVSplitter->addWidget(fHelpArea);
#endif
  
  
#if QT_VERSION >= 0x040000
  vLayout->addWidget(helpWidget);
  vLayout->addWidget(fHelpVSplitter,1);
#endif
  
  fHelpTBWidget->setMinimumSize(50,50);
  fHelpTBWidget->setSizePolicy (QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));
  // set the splitter size
#if QT_VERSION >= 0x040000
  QList<int> list;
#else
  QValueList<int> list;
#endif
  list.append( 50 );
  list.append( 50 );
  fHelpVSplitter->setSizes(list);
  
#if QT_VERSION >= 0x040000
  helpWidget->setLayout(helpLayout);
  fHelpTBWidget->setLayout(vLayout);
#endif  
}


/** Create the Cout ToolBox Widget
 */
void G4UIQt::CreateCoutTBWidget(
) 
{
#if QT_VERSION >= 0x040000
  QVBoxLayout *layoutCoutTB = new QVBoxLayout();
#else
  QVBoxLayout *layoutCoutTB = new QVBoxLayout(fCoutTBWidget);
#endif

  fCoutTBTextArea = new QTextEdit(fCoutTBWidget);
  fCoutFilter = new QLineEdit(fCoutTBWidget);
  QLabel* coutFilterLabel = new QLabel("Filter : ",fCoutTBWidget);

  QPushButton *coutTBClearButton = new QPushButton("clear",fCoutTBWidget);
  connect(coutTBClearButton, SIGNAL(clicked()), SLOT(ClearButtonCallback()));
  connect(fCoutFilter, SIGNAL(textEdited ( const QString &)), SLOT(CoutFilterCallback( const QString &)));

  fCoutTBTextArea->setReadOnly(true);

  QWidget* coutButtonWidget = new QWidget(fCoutTBWidget);
  QHBoxLayout* layoutCoutTBButtons = new QHBoxLayout(coutButtonWidget);
  layoutCoutTBButtons->addWidget(coutTBClearButton);
  layoutCoutTBButtons->addWidget(coutFilterLabel);
  layoutCoutTBButtons->addWidget(fCoutFilter);

  layoutCoutTB->addWidget(fCoutTBTextArea);
  layoutCoutTB->addWidget(coutButtonWidget);

#if QT_VERSION >= 0x040000
  fCoutTBWidget->setLayout(layoutCoutTB);
#endif
}


/** Create the VisParameters ToolBox Widget
 */
void G4UIQt::CreateVisParametersTBWidget(
) 
{
}


/** Create the ViewComponents ToolBox Widget
 */
void G4UIQt::CreateViewComponentsTBWidget(
) 
{
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
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::AddTabWidget %d %d\n",sizeX, sizeY);
#endif

  if (fTabWidget == NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIQt::AddTabWidget +++++\n");
#endif
    fTabWidget = new G4QTabWidget(fMyVSplitter);
#if QT_VERSION >= 0x040500
    fTabWidget->setTabsClosable (true); 
#endif
    
#if QT_VERSION >= 0x040200
    fTabWidget->setUsesScrollButtons (true);
#endif
    
    fTabWidget->setSizePolicy (QSizePolicy(QSizePolicy::Maximum,QSizePolicy::Maximum));
    
    QSizePolicy policy = fTabWidget->sizePolicy();
#if QT_VERSION < 0x040000
    policy.setHorStretch(1);
    policy.setVerStretch(1);
#else
    policy.setHorizontalStretch(1);
    policy.setVerticalStretch(1);
#endif
    fTabWidget->setSizePolicy(policy);
    
#if QT_VERSION >= 0x040500
    connect(fTabWidget,   SIGNAL(tabCloseRequested(int)), this, SLOT(TabCloseCallback(int)));
#endif
    connect(fTabWidget, SIGNAL(currentChanged ( int ) ), SLOT(UpdateTabWidget(int))); 
  }

  fLastQTabSizeX = sizeX;
  fLastQTabSizeY = sizeY;

  if (!aWidget) {
    return false;
  }

  // Remove QLabel 

  // L.Garnier 26/05/2010 : not exactly the same in qt3. Could cause some
  // troubles
#if QT_VERSION >= 0x040000
  if ( fMyVSplitter->indexOf(fEmptyViewerTabLabel) != -1) {
#endif

#if QT_VERSION < 0x040000
    fEmptyViewerTabLabel->reparent(0,0,QPoint(0,0));  
    fEmptyViewerTabLabel->hide();
    delete fEmptyViewerTabLabel;
    fEmptyViewerTabLabel = NULL;

#else
    fEmptyViewerTabLabel->hide();
    fEmptyViewerTabLabel->setParent(NULL);
    delete fEmptyViewerTabLabel;
    fEmptyViewerTabLabel = NULL;

    fMyVSplitter->addWidget(fTabWidget);
#endif

#if QT_VERSION < 0x040000
    aWidget->reparent(fTabWidget,0,QPoint(0,0));  
#else
    aWidget->setParent(fTabWidget);
#endif
#if QT_VERSION >= 0x040000
  }
#endif


#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::AddTabWidget ADD %d %d + %d %d---------------------------------------------------\n",sizeX, sizeY,sizeX-fTabWidget->width(),sizeY-fTabWidget->height());
#endif

  if (fMainWindow->isVisible()) {

    // get the size of the tabbar
    int tabBarX = 0;
    int tabBarY = 0;
    if (fTabWidget->count() >0) {
#if QT_VERSION < 0x040000
      tabBarX = fTabWidget->width()-fTabWidget->page(0)->width();
      tabBarY = fTabWidget->height()-fTabWidget->page(0)->height();
#else
      tabBarX = fTabWidget->width()-fTabWidget->widget(0)->width();
      tabBarY = fTabWidget->height()-fTabWidget->widget(0)->height();
#endif
    }

    fMainWindow->resize(tabBarX+fMainWindow->width()+sizeX-fTabWidget->width(),tabBarY+fMainWindow->height()+sizeY-fTabWidget->height());
  }

  // Problems with resize. The widgets are not realy drawn at this step,
  // then we have to force them on order to check the size

#if QT_VERSION < 0x040000
  fTabWidget->insertTab(aWidget,name,fTabWidget->count());
#else
  fTabWidget->insertTab(fTabWidget->count(),aWidget,name);
#endif
  
#if QT_VERSION < 0x040000
  fTabWidget->setCurrentPage(fTabWidget->count()-1);
#else
  fTabWidget->setCurrentIndex(fTabWidget->count()-1);
#endif

  // Set visible
#if QT_VERSION >= 0x040000
 #if QT_VERSION >= 0x040200
   fTabWidget->setLastTabCreated(fTabWidget->currentIndex());
 #else
   fTabWidget->setLastTabCreated(fTabWidget->currentIndex());
 #endif
#else
  fTabWidget->setLastTabCreated(fTabWidget->currentPageIndex());
#endif
  
  return true;
}


void G4UIQt::UpdateTabWidget(int tabNumber) {
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::UpdateTabWidget %d\n",tabNumber);
#endif
  if (  fTabWidget == NULL) {
    fTabWidget = new G4QTabWidget;
  }
  

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::UpdateTabWidget CALL REPAINT tabGL\n");
#endif

#if QT_VERSION < 0x040000
  fTabWidget->setCurrentPage(tabNumber);
#else
  fTabWidget->setCurrentIndex(tabNumber);
#endif

  // Send this signal to unblock graphic updates !
  fTabWidget->setTabSelected(false);

#if QT_VERSION >= 0x040000
 #if QT_VERSION >= 0x040200
  fTabWidget->setVisible(true);
 #else
  fTabWidget->show();
 #endif
#else
  fTabWidget->show();
#endif

  // This will send a paintEvent to OGL Viewers
  fTabWidget->setTabSelected(true);

#if QT_VERSION < 0x040000
  QApplication::sendPostedEvents () ;
#else
  QCoreApplication::sendPostedEvents () ;
#endif

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::UpdateTabWidget END\n");
#endif
}


/** Send resize event to all tabs
 */
void G4UIQt::ResizeTabWidget( QResizeEvent* e) {
  for (G4int a=0;a<fTabWidget->count() ;a++) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIQt::ResizeTabWidget +++++++++++++++++++++++++++++++++++++++\n");
#endif
#if QT_VERSION < 0x040000
    fTabWidget->page(a)->resize(e->size());
#else
    fTabWidget->widget(a)->resize(e->size());
#endif
  }
}


/**   Start the Qt main loop
*/
G4UIsession* G4UIQt::SessionStart (
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt SessionStart\n");
#endif

  G4Qt* interactorManager = G4Qt::getInstance ();

  Prompt("Session :");
  exitSession = false;

  if (fEmptyViewerTabLabel != NULL) {
    bool visible = false;
    if (fTabWidget != NULL) {
      if (fTabWidget->isVisible()) {
        visible = true;
      }
    }
  }

#if QT_VERSION >= 0x040000
  #if QT_VERSION >= 0x040200
      fMainWindow->setVisible(true);
  #else
      fMainWindow->show();
  #endif
#else
      fMainWindow->show();
#endif
  // get the size of the tabbar
  int tabBarX = 0;
  int tabBarY = 0;

  if (fTabWidget != NULL) {
#if QT_VERSION < 0x040000
    tabBarX = -fTabWidget->page(0)->width();
    tabBarY = -fTabWidget->page(0)->height();
#else
    tabBarX = -fTabWidget->widget(0)->width();
    tabBarY = -fTabWidget->widget(0)->height();
#endif
  }
  fMainWindow->resize(tabBarX+fMainWindow->width()+fLastQTabSizeX,tabBarY+fMainWindow->height()+fLastQTabSizeY);

#if QT_VERSION < 0x040000
  QApplication::sendPostedEvents () ;
#else
  QCoreApplication::sendPostedEvents () ;
#endif

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt SessionStart2\n");
#endif
  interactorManager->DisableSecondaryLoop (); // TO KEEP
  if ((QApplication*)interactorManager->GetMainInteractor())
    ((QApplication*)interactorManager->GetMainInteractor())->exec();

  // on ne passe pas le dessous ? FIXME ????
  // je ne pense pas 13/06

  //   void* event; // TO KEEP
  //   while((event = interactorManager->GetEvent())!=NULL) {  // TO KEEP
  //     interactorManager->DispatchEvent(event); // TO KEEP
  //     if(exitSession==true) break; // TO KEEP
  //   } // TO KEEP

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
 G4String aState
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
  void* event; // TO KEEP
  while((event = interactorManager->GetEvent())!=NULL) {  // TO KEEP
    interactorManager->DispatchEvent(event); // TO KEEP
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
 G4String aString
 )
{
  if (!aString) return 0;
  
  QStringList newStr;
  
  // Add to stringList
#if QT_VERSION < 0x040000
  newStr = QStringList(QString((char*)aString.data()).simplifyWhiteSpace());
#else
  newStr = QStringList(QString((char*)aString.data()).trimmed());
#endif
  fG4cout += newStr;
 
#if QT_VERSION >= 0x040000
  QStringList result = newStr.filter(fCoutFilter->text());
#else
  //L. Garnier : in qt3 filter will does nothing
  QStringList result = "";
#endif

  if (result.join("\n").isEmpty()) {
    return 0;
  }
  fCoutTBTextArea->append(result.join("\n"));
  fCoutTBTextArea->repaint();

#if QT_VERSION < 0x040000
  fCoutTBTextArea->verticalScrollBar()->setValue(fCoutTBTextArea->verticalScrollBar()->maxValue());
#else
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
#endif

  return 0;
}


/**
   Receive a cerr from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIQt::ReceiveG4cerr (
 G4String aString
)
{
  if (!aString) return 0;

  QStringList newStr;

  // Add to stringList
#if QT_VERSION < 0x040000
  newStr = QStringList(QString((char*)aString.data()).simplifyWhiteSpace());
#else
  newStr = QStringList(QString((char*)aString.data()).trimmed());
#endif
  fG4cout += newStr;
 
#if QT_VERSION < 0x040000
  //L. Garnier : in qt3 filter will does nothing
  QStringList result = "";
#else
  QStringList result = newStr.filter(fCoutFilter->text());
#endif

#if QT_VERSION < 0x040000
  QColor previousColor = fCoutTBTextArea->color();
  fCoutTBTextArea->setColor(Qt::red);
  fCoutTBTextArea->append(result.join("\n"));
  fCoutTBTextArea->setColor(previousColor);
  fCoutTBTextArea->verticalScrollBar()->setValue(fCoutTBTextArea->verticalScrollBar()->maxValue());
#else
  QColor previousColor = fCoutTBTextArea->textColor();
  fCoutTBTextArea->setTextColor(Qt::red);
  fCoutTBTextArea->append(result.join("\n"));
  fCoutTBTextArea->setTextColor(previousColor);
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
#endif
  fCoutTBTextArea->repaint();
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

#if QT_VERSION < 0x040000
  QPopupMenu *fileMenu = new QPopupMenu( fMainWindow);
  fMainWindow->menuBar()->insertItem( aLabel, fileMenu );
#else
  QMenu *fileMenu = new QMenu(aLabel);
  fMainWindow->menuBar()->insertMenu(fMainWindow->menuBar()->actions().last(),fileMenu); 
#endif

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

#if QT_VERSION < 0x040000
  QPopupMenu *parent = (QPopupMenu*)GetInteractor(aMenu);
#else
  QMenu *parent = (QMenu*)GetInteractor(aMenu);
#endif

  if(parent==NULL) return;
  
  QSignalMapper *signalMapper = new QSignalMapper(this);
#if QT_VERSION < 0x030200
  QAction *action = new QAction(QString(aLabel),QString(aLabel),QKeySequence(),signalMapper, SLOT(map()));
  action->addTo(parent);
 connect(action,SIGNAL(activated()),signalMapper,SLOT(map()));

#elif QT_VERSION < 0x040000
  QAction *action = new QAction(QString(aLabel),QKeySequence(),signalMapper, SLOT(map()));
  action->addTo(parent);
 connect(action,SIGNAL(activated()),signalMapper,SLOT(map()));

#else
  QAction *action = parent->addAction(aLabel, signalMapper, SLOT(map()));

#endif
  connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ButtonCallback(const QString&)));
  signalMapper->setMapping(action, QString(aCommand));
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
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::ActivateCommand found : %s \n",targetCom.data());
#endif
  if (targetCom != "") {
    OpenHelpTreeOnCommand(targetCom.data());
  }

#if QT_VERSION < 0x040000
  fToolBox->setCurrentItem(fHelpTBWidget);
#else
  fToolBox->setCurrentWidget(fHelpTBWidget);
#endif
}



/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIQt::InitHelpTree()
{

  if (! fHelpTreeWidget ) {
#if QT_VERSION < 0x040000
    fHelpTreeWidget = new QListView(fHelpVSplitter);
#else
    fHelpTreeWidget = new QTreeWidget();
#endif
  }


  // build widget
#if QT_VERSION < 0x040000
  fHelpTreeWidget->setSelectionMode(QListView::Single);
  fHelpTreeWidget->setRootIsDecorated(true);
  fHelpTreeWidget->addColumn("Command");
  fHelpTreeWidget->header()->setResizeEnabled(FALSE,1);
#else
  fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  QStringList labels;
  labels << QString("Command");
  fHelpTreeWidget->setHeaderLabels(labels);
#endif


#if QT_VERSION < 0x040000
  connect(fHelpTreeWidget, SIGNAL(selectionChanged ()),this, SLOT(HelpTreeClicCallback()));  
  connect(fHelpTreeWidget, SIGNAL(doubleClicked (QListViewItem*)),this, SLOT(HelpTreeDoubleClicCallback()));
#else
  connect(fHelpTreeWidget, SIGNAL(itemSelectionChanged ()),this, SLOT(HelpTreeClicCallback()));  
  connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked (QTreeWidgetItem*,int)),this, SLOT(HelpTreeDoubleClicCallback()));  
#endif

}
/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIQt::FillHelpTree()
{
  if (! fHelpTreeWidget ) {
    InitHelpTree();
  }

  QString searchText = fHelpLine->text();

  if (searchText =="") {
    // clear old help tree
    //    fHelpTreeWidget->clear();
#if QT_VERSION < 0x040000
    fHelpTreeWidget->removeColumn(1);
    fHelpTreeWidget->removeColumn(0);
#endif
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
#if QT_VERSION < 0x040000
  QListViewItem * newItem = NULL;
#else
  QTreeWidgetItem * newItem = NULL;
#endif
  QString commandText = "";
  for (int a=0;a<treeSize;a++) {
    // Creating new item
    newItem = NULL;

#if QT_VERSION < 0x040000
    commandText = QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).simplifyWhiteSpace();
#else
    commandText = QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).trimmed();
#endif

    // if already exist, don't create it !
#if QT_VERSION < 0x040000
    QListViewItem* tmpAddItem = fHelpTreeWidget->firstChild();
    while (tmpAddItem != 0) {
      if (!newItem) {
        newItem = FindTreeItem(tmpAddItem,commandText);
      }
      tmpAddItem = tmpAddItem->nextSibling();
    }
#else
    for (int b=0;b<fHelpTreeWidget->topLevelItemCount();b++) {
      if (!newItem)
        newItem = FindTreeItem(fHelpTreeWidget->topLevelItem(b),commandText);
    }
#endif

    if (newItem == NULL) {
      
#if QT_VERSION < 0x040000
      newItem = new QListViewItem(fHelpTreeWidget);
#else
      newItem = new QTreeWidgetItem(fHelpTreeWidget);
#endif
      newItem->setText(0,GetShortCommandPath(commandText));
    }

    // look for childs
    CreateChildTree(newItem,treeTop->GetTree(a+1));
  }

}



/**   Fill the Help Tree Widget
   @param aParent : parent item to fill
   @param aCommandTree : commandTree node associate with this part of the Tree
*/
#if QT_VERSION < 0x040000
void G4UIQt::CreateChildTree(
 QListViewItem *aParent
,G4UIcommandTree *aCommandTree
#else
void G4UIQt::CreateChildTree(
 QTreeWidgetItem *aParent
,G4UIcommandTree *aCommandTree
#endif
)
{
  if (aParent == NULL) return;
  if (aCommandTree == NULL) return;


  // Creating new item
#if QT_VERSION < 0x040000
  QListViewItem * newItem;
#else
  QTreeWidgetItem * newItem;
#endif

  QString commandText = "";
  // Get the Sub directories
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {

#if QT_VERSION < 0x040000
    commandText = QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).simplifyWhiteSpace();
#else
    commandText = QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).trimmed();
#endif
    
    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
#if QT_VERSION < 0x040000
      newItem = new QListViewItem(aParent);
#else
      newItem = new QTreeWidgetItem(aParent);
#endif
      newItem->setText(0,GetShortCommandPath(commandText));
    }
    CreateChildTree(newItem,aCommandTree->GetTree(a+1));
  }



  // Get the Commands

  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    QStringList stringList;
#if QT_VERSION < 0x040000
    commandText = QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).simplifyWhiteSpace();
#else
    commandText = QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed();
#endif

    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
#if QT_VERSION < 0x040000
      newItem = new QListViewItem(aParent);
      newItem->setText(0,GetShortCommandPath(commandText));
      newItem->setOpen(false);
      
#else
      newItem = new QTreeWidgetItem(aParent);
      newItem->setText(0,GetShortCommandPath(commandText));
#if QT_VERSION < 0x040202
      fHelpTreeWidget->setItemExpanded(newItem,false); 
#else
      newItem->setExpanded(false);
#endif
#endif
    }
  }
}

 
/** Find a treeItemWidget in the help tree
    @param aCommand item's String to look for
    @return item if found, NULL if not
*/
#if QT_VERSION < 0x040000
QListViewItem* G4UIQt::FindTreeItem(
 QListViewItem *aParent
#else
QTreeWidgetItem* G4UIQt::FindTreeItem(
 QTreeWidgetItem *aParent
#endif
,const QString& aCommand
)
{
  if (aParent == NULL) return NULL;

  // Suppress last "/"
  QString myCommand = aCommand;
  
#if QT_VERSION < 0x040000
  if (myCommand.findRev("/") == ((int)myCommand.length()-1)) {
    myCommand = myCommand.left(myCommand.length()-1);
#else
  if (myCommand.lastIndexOf("/") == (myCommand.size()-1)) {
    myCommand = myCommand.left(myCommand.size()-1);
#endif
  }

  if (GetLongCommandPath(aParent) == myCommand)
    return aParent;
  
#if QT_VERSION < 0x040000
  QListViewItem * tmp = NULL;
  QListViewItem* tmpItem = aParent->firstChild();
    while (tmpItem != 0) {
      if (!tmp)
        tmp = FindTreeItem(tmpItem,myCommand);
      tmpItem = tmpItem->nextSibling();
    }
#else
  QTreeWidgetItem * tmp = NULL;
  for (int a=0;a<aParent->childCount();a++) {
    if (!tmp)
      tmp = FindTreeItem(aParent->child(a),myCommand);
  }
#endif
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



/**  Implement G4VBasicShell vurtual function
 */
G4bool G4UIQt::GetHelpChoice(
 G4int&
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::GetHelpChoice SHOULD NEVER GO HERE");
#endif
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
#if QT_VERSION < 0x040000
        // count rows...
        QListViewItem* tmpItem = fHistoryTBTableList->firstChild();
        int selection = -1;
        int index = 0;
        while (tmpItem != 0) {
          if (tmpItem == fHistoryTBTableList->selectedItem()) {
            selection = index;
          }
          index ++;
          tmpItem = tmpItem->nextSibling();
        }
        if (fHistoryTBTableList->childCount()) {
          if (selection == -1) {
            selection = fHistoryTBTableList->childCount()-1;
          } else {
            if (e->key() == (Qt::Key_Down)) {
              if (selection <(fHistoryTBTableList->childCount()-1))
                selection++;
            } else if (e->key() == (Qt::Key_PageDown)) {
              selection = fHistoryTBTableList->childCount()-1;
#else
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
#endif
            } else if (e->key() == (Qt::Key_Up)) {
              if (selection >0)
                selection --;
            } else if (e->key() == (Qt::Key_PageUp)) {
              selection = 0;
            }
          }
          fHistoryTBTableList->clearSelection();
#if QT_VERSION < 0x040000
          QListViewItem* tmpItem = fHistoryTBTableList->firstChild();
          int index = 0;
          while (tmpItem != 0) {
            if (index == selection) {
              tmpItem->setSelected(true);
              fHistoryTBTableList->setCurrentItem(tmpItem);
          }
          index ++;
          tmpItem = tmpItem->nextSibling();
        }
#else
#if QT_VERSION < 0x040202
          fHistoryTBTableList->setItemSelected(fHistoryTBTableList->item(selection),true);
#else
          fHistoryTBTableList->item(selection)->setSelected(true);
#endif      
          fHistoryTBTableList->setCurrentItem(fHistoryTBTableList->item(selection));
#endif
        }
        moveCommandCursor = true;
      } else if (e->key() == (Qt::Key_Tab)) {
#if QT_VERSION < 0x040000
        G4String ss = Complete(fCommandArea->text().ascii());
#else
        G4String ss = Complete(fCommandArea->text().toStdString().c_str());
#endif
        fCommandArea->setText((char*)(ss.data()));

        // do not pass by parent, it will disable widget tab focus !
        return true;
#if QT_VERSION >= 0x040000
        // L.Garnier : MetaModifier is CTRL for MAC, but I don't want to put a MAC 
        // specific #ifdef
      } else if (((e->modifiers () == Qt::ControlModifier) || (e->modifiers () == Qt::MetaModifier)) && (e->key() == Qt::Key_A)) {
       fCommandArea->home(false);
       return true;
      } else if (((e->modifiers () == Qt::ControlModifier) || (e->modifiers () == Qt::MetaModifier)) && (e->key() == Qt::Key_E)) {
       fCommandArea->end(false);
       return true;
#endif
      }
    }
  }
  bool res= false;
  // change cursor position if needed
  if (moveCommandCursor == true) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIQt::eventFilter setCursor Position\n");
#endif
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
)
{
}

/**   Callback call when "click on a menu entry.<br>
   Send the associated command to geant4
*/
void G4UIQt::CommandEnteredCallback (
)
{
#if QT_VERSION < 0x040000
  G4String command (fCommandArea->text().ascii());
  if (fCommandArea->text().simplifyWhiteSpace() != "") {

    QListViewItem *newItem = new QListViewItem(fHistoryTBTableList);
    newItem->setText(0,fCommandArea->text());
    fHistoryTBTableList->insertItem(newItem);
    // now we have to arrange 
    QListViewItem *temp= fHistoryTBTableList->lastItem();
    for (int i=0; i<fHistoryTBTableList->childCount()-1;i++) {
      fHistoryTBTableList->takeItem(temp);
      fHistoryTBTableList->insertItem(temp);
      temp= fHistoryTBTableList->lastItem();
    }
#else
  G4String command (fCommandArea->text().toStdString().c_str());
  if (fCommandArea->text().trimmed() != "") {
    fHistoryTBTableList->addItem(fCommandArea->text());
#endif
    fHistoryTBTableList->clearSelection();
    fHistoryTBTableList->setCurrentItem(NULL);
    fCommandArea->setText("");

    G4Qt* interactorManager = G4Qt::getInstance ();
    if (interactorManager) { 
      interactorManager->FlushAndWaitExecution();
    }
    if (command(0,4) != "help") {
      ApplyShellCommand (command,exitSession,exitPause);
    } else {
      ActivateCommand(command);
    }
    // Rebuild help tree
    FillHelpTree();

    if(exitSession==true) 
      SessionTerminate();
  }
}


/**   Callback call when "enter" clicked on the command zone.<br>
   Send the command to geant4
   @param aCommand
*/
void G4UIQt::ButtonCallback (
 const QString& aCommand
)
{
#if QT_VERSION < 0x040000
  G4String ss = G4String(aCommand.ascii());
#else
  G4String ss = G4String(aCommand.toStdString().c_str());
#endif
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
#if QT_VERSION < 0x040000
  QListViewItem* item =  NULL;
#else
  QTreeWidgetItem* item =  NULL;
#endif
  if (!fHelpTreeWidget)
    return ;

  if (!fHelpArea)
    return;
  
#if QT_VERSION < 0x040000
  item =fHelpTreeWidget->selectedItem();
#else
  QList<QTreeWidgetItem *> list =fHelpTreeWidget->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
#endif
  if (!item)
    return;
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();

#if QT_VERSION < 0x040000
  std::string itemText = GetLongCommandPath(item).ascii();
#else
  std::string itemText = GetLongCommandPath(item).toStdString();
#endif
  
  G4UIcommand* command = treeTop->FindPath(itemText.c_str());

  if (command) {
#if QT_VERSION >= 0x040000
#if QT_VERSION < 0x040200
    fHelpArea->clear();
    fHelpArea->append(GetCommandList(command));
#else
    fHelpArea->setText(GetCommandList(command));
#endif
#else
    fHelpArea->setText(GetCommandList(command));
#endif
  } else {  // this is a command
    G4UIcommandTree* path = treeTop->FindCommandTree(itemText.c_str());
    if ( path) {
      // this is not a command, this is a sub directory
      // We display the Title
#if QT_VERSION >= 0x040000
#if QT_VERSION < 0x040200
      fHelpArea->clear();
      fHelpArea->append(path->GetTitle().data());
#else
      fHelpArea->setText(path->GetTitle().data());
#endif
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

#if QT_VERSION < 0x040000
  QListViewItem* item =  NULL;
#else
  QTreeWidgetItem* item =  NULL;
#endif
  if (!fHelpTreeWidget)
    return ;

  if (!fHelpArea)
    return;
  
#if QT_VERSION < 0x040000
  item =fHelpTreeWidget->selectedItem();
#else
  QList<QTreeWidgetItem *> list =fHelpTreeWidget->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
#endif
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
#if QT_VERSION < 0x040000
  QListViewItem* item =  NULL;
#else
  QListWidgetItem* item =  NULL;
#endif
  if (!fHistoryTBTableList)
    return ;

  
#if QT_VERSION < 0x040000
  item =fHistoryTBTableList->selectedItem();
#else
  QList<QListWidgetItem *> list =fHistoryTBTableList->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
#endif
  if (!item)
    return;
#if QT_VERSION < 0x040000
  fCommandArea->setText(item->text(0));
#else
  fCommandArea->setText(item->text());
#endif
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIQt::CommandHistoryCallback change text\n");
#endif
}


void G4UIQt::CoutFilterCallback(
#if QT_VERSION < 0x040000
const QString &) {
#else
const QString & text) {
#endif

#if QT_VERSION < 0x040000
  QStringList result = "";
#else
  QStringList result = fG4cout.filter(text);
  fCoutTBTextArea->setPlainText(result.join("\n"));
#endif

  fCoutTBTextArea->repaint();
#if QT_VERSION < 0x040000
  fCoutTBTextArea->verticalScrollBar()->setValue(fCoutTBTextArea->verticalScrollBar()->maxValue());
#else
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
#endif

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
#if QT_VERSION < 0x040000
    fHelpTreeWidget->removeColumn(1);
    fHelpTreeWidget->removeColumn(0);
#endif

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
#if QT_VERSION < 0x040000
  fHelpTreeWidget->removeColumn(1);
  fHelpTreeWidget->removeColumn(0);
#endif

  // look for new items

  int tmp = 0;
#if QT_VERSION < 0x040000
  int multFactor = 1000; // factor special for having double keys in Qt3
  int doubleKeyAdd = 0;  // decay for doubleKeys in Qt3
#endif

  QMap<int,QString> commandResultMap;
  QMap<int,QString> commandChildResultMap;

  for (int a=0;a<treeSize;a++) {
    G4UIcommand* command = treeTop->FindPath(treeTop->GetTree(a+1)->GetPathName().data());
#if QT_VERSION > 0x040000
    tmp = GetCommandList (command).count(searchText,Qt::CaseInsensitive);
#else
    tmp = GetCommandList (command).contains(searchText,false);
#endif
    if (tmp >0) {
#if QT_VERSION > 0x040000
      commandResultMap.insertMulti(tmp,QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()));
#else // tricky thing for Qt3...
      doubleKeyAdd = 0;
      while (commandResultMap.find( tmp*multFactor+doubleKeyAdd) != commandResultMap.end()) {
        doubleKeyAdd ++;
      }
      commandResultMap.insert( tmp*multFactor+doubleKeyAdd,QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()) );
#endif
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(treeTop->GetTree(a+1),searchText);
    // insert new childs
    if (!commandChildResultMap.empty()) {
#if QT_VERSION > 0x040000
      QMap<int,QString>::const_iterator i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insertMulti(i.key(),i.value());
#else // tricky thing for Qt3...
      QMap<int,QString>::const_iterator i = commandChildResultMap.begin();
      while (i != commandChildResultMap.end()) {
        doubleKeyAdd = 0;
        while (commandResultMap.find( i.key()*multFactor+doubleKeyAdd) != commandResultMap.end()) {
          doubleKeyAdd ++;
        }
        commandResultMap.insert(i.key()*multFactor+doubleKeyAdd,i.data());
#endif
        i++;
      }
      commandChildResultMap.clear();
    }
  }

  // build new help tree
#if QT_VERSION < 0x040000
  fHelpTreeWidget->setSelectionMode(QListView::Single);
  fHelpTreeWidget->setRootIsDecorated(true);
  fHelpTreeWidget->addColumn("Command");
  fHelpTreeWidget->addColumn("Match");
  //  fHelpTreeWidget->header()->setResizeEnabled(FALSE,1);
#else
  fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  fHelpTreeWidget->setColumnCount(2);
  QStringList labels;
  labels << QString("Command") << QString("Match");
  fHelpTreeWidget->setHeaderLabels(labels);
#endif

  if (commandResultMap.empty()) {
#if QT_VERSION < 0x040200
    fHelpArea->clear();
    fHelpArea->append("No match found");
#else
    fHelpArea->setText("No match found");
#endif
    return;
  }

#if QT_VERSION > 0x040000
  QMap<int,QString>::const_iterator i = commandResultMap.constEnd();
#else
  QMap<int,QString>::const_iterator i = commandResultMap.end();
#endif
  i--;
  // 10 maximum progress values
  float multValue = 10.0/(float)(i.key());
  QString progressChar = "|";
  QString progressStr = "|";

#if QT_VERSION < 0x040000
  QListViewItem * newItem;
#else
  QTreeWidgetItem * newItem;
#endif
  bool end = false;
  while (!end) {
#if QT_VERSION > 0x040000
    if (i == commandResultMap.constBegin()) {
#else
    if (i == commandResultMap.begin()) {
#endif
      end = true;
    }
    for(int a=0;a<int(i.key()*multValue);a++) {
      progressStr += progressChar;
    }
#if QT_VERSION < 0x040000
    newItem = new QListViewItem(fHelpTreeWidget);
    QString commandStr = i.data().simplifyWhiteSpace();
#else
    newItem = new QTreeWidgetItem(fHelpTreeWidget);
    QString commandStr = i.value().trimmed();
#endif

#if QT_VERSION < 0x040000
    if (commandPath.find("/") == 0) {
      commandStr = commandStr.right(commandStr.length()-1);
#else
    if (commandStr.indexOf("/") == 0) {
      commandStr = commandStr.right(commandStr.size()-1);
#endif
    }
      
    newItem->setText(0,commandStr);
    newItem->setText(1,progressStr);
    
#if QT_VERSION >= 0x040200
    newItem->setForeground ( 1, QBrush(Qt::blue) );
#endif
    progressStr = "|";
    i--;
  }
  // FIXME :  to be checked on Qt3
#if QT_VERSION < 0x040000
  fHelpTreeWidget->setColumnWidthMode (1,QListView::Maximum);
  fHelpTreeWidget->setSorting(1,false);
#else
  fHelpTreeWidget->resizeColumnToContents (0);
  fHelpTreeWidget->sortItems(1,Qt::DescendingOrder);
  //  fHelpTreeWidget->setColumnWidth(1,10);//resizeColumnToContents (1);
#endif
}




QMap<int,QString> G4UIQt::LookForHelpStringInChildTree(
 G4UIcommandTree *aCommandTree
,const QString & text
 )
{
  QMap<int,QString> commandResultMap;
  if (aCommandTree == NULL) return commandResultMap;
  
#if QT_VERSION < 0x040000
  int multFactor = 1000; // factor special for having double keys in Qt3
  int doubleKeyAdd = 0;  // decay for doubleKeys in Qt3
#endif

  // Get the Sub directories
  int tmp = 0;
  QMap<int,QString> commandChildResultMap;
  
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {
    const G4UIcommand* command = aCommandTree->GetGuidance();
#if QT_VERSION > 0x040000
    tmp = GetCommandList (command).count(text,Qt::CaseInsensitive);
#else
    tmp = GetCommandList (command).contains(text,false);
#endif
    if (tmp >0) {
#if QT_VERSION > 0x040000
      commandResultMap.insertMulti(tmp,QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()));
#else // tricky thing for Qt3...
      doubleKeyAdd = 0;
      while (commandResultMap.find( tmp*multFactor+doubleKeyAdd) != commandResultMap.end()) {
        doubleKeyAdd ++;
      }
      commandResultMap.insert(tmp*multFactor+doubleKeyAdd,QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()));
#endif
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(aCommandTree->GetTree(a+1),text);
    
    if (!commandChildResultMap.empty()) {
      // insert new childs
#if QT_VERSION > 0x040000
      QMap<int,QString>::const_iterator i = commandChildResultMap.constBegin();
      while (i != commandChildResultMap.constEnd()) {
        commandResultMap.insertMulti(i.key(),i.value());
#else // tricky thing for Qt3...
      QMap<int,QString>::const_iterator i = commandChildResultMap.begin();
      while (i != commandChildResultMap.end()) {
        doubleKeyAdd = 0;
        while (commandResultMap.find( i.key()*multFactor+doubleKeyAdd) != commandResultMap.end()) {
          doubleKeyAdd ++;
        }
        commandResultMap.insert(i.key()*multFactor+doubleKeyAdd,i.data());
#endif
        i++;
      }
      commandChildResultMap.clear();
    }
  }
  // Get the Commands
  
  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    const G4UIcommand* command = aCommandTree->GetCommand(a+1);
#if QT_VERSION > 0x040000
    tmp = GetCommandList (command).count(text,Qt::CaseInsensitive);
#else
    tmp = GetCommandList (command).contains(text,false);
#endif
    if (tmp >0) {
#if QT_VERSION > 0x040000
      commandResultMap.insertMulti(tmp,QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()));
#else // tricky thing for Qt3...
      doubleKeyAdd = 0;
      while (commandResultMap.find( tmp*multFactor+doubleKeyAdd) != commandResultMap.end()) {
        doubleKeyAdd ++;
      }
      commandResultMap.insert(tmp*multFactor+doubleKeyAdd,QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()));
#endif
    }
    
  }
  return commandResultMap;
}

  
QString G4UIQt::GetShortCommandPath(
QString commandPath
)
{
#if QT_VERSION < 0x040000
  if (commandPath.find("/") == 0) {
    commandPath = commandPath.right(commandPath.length()-1);
#else
  if (commandPath.indexOf("/") == 0) {
    commandPath = commandPath.right(commandPath.size()-1);
#endif
  }

#if QT_VERSION < 0x040000
  commandPath = commandPath.right(commandPath.length()-commandPath.findRev("/",-2)-1);
#else
  commandPath = commandPath.right(commandPath.size()-commandPath.lastIndexOf("/",-2)-1);
#endif
 
#if QT_VERSION < 0x040000
  if (commandPath.findRev("/") == ((int)commandPath.length()-1)) {
    commandPath = commandPath.left(commandPath.length()-1);
 }
#else
 if (commandPath.lastIndexOf("/") == (commandPath.size()-1)) {
    commandPath = commandPath.left(commandPath.size()-1);
 }
#endif

 return commandPath;
}


QString G4UIQt::GetLongCommandPath(
#if QT_VERSION < 0x040000
 QListViewItem* item
#else
 QTreeWidgetItem* item
#endif
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

G4QTabWidget::G4QTabWidget(
QSplitter*& split
):QTabWidget(split)
 ,tabSelected(false)
 ,lastCreated(-1)
{
}

G4QTabWidget::G4QTabWidget(
):QTabWidget()
 ,tabSelected(false)
 ,lastCreated(-1)
{
}


  
#if QT_VERSION >= 0x040500
void G4UIQt::TabCloseCallback(int a){
#else
void G4UIQt::TabCloseCallback(int){
#endif
#if QT_VERSION >= 0x040500
  QWidget* temp = fTabWidget->widget(a);
  fTabWidget->removeTab (a);

  delete temp;

  if (fTabWidget->count() == 0) {
    if (fEmptyViewerTabLabel == NULL) {
#if QT_VERSION < 0x040000
      fEmptyViewerTabLabel = new QLabel(fToolBox,"         If you want to have a Viewer, please use /vis/open commands. ");
#else
      fEmptyViewerTabLabel = new QLabel("         If you want to have a Viewer, please use /vis/open commands. ");
#endif
    }

    fMyVSplitter->addWidget(fEmptyViewerTabLabel);
    fMyVSplitter->show();
    fEmptyViewerTabLabel->show();
    fTabWidget->setParent(0);
#if QT_VERSION >= 0x040000
  #if QT_VERSION >= 0x040200
      fTabWidget->setVisible(false);
  #else
      fMainWindow->hide();
  #endif
#else
      fMainWindow->hide();
#endif
    delete fTabWidget;
    fTabWidget = NULL;
  }
#endif
}


void G4UIQt::ToolBoxActivated(int a){

#if QT_VERSION < 0x040000
  if (fToolBox->item(a) == fHelpTBWidget) {
#else
  if (fToolBox->widget(a) == fHelpTBWidget) {
#endif
    // Rebuild the help tree
    FillHelpTree();
  }
}

void G4QTabWidget::paintEvent(
QPaintEvent *
)
{

#if QT_VERSION < 0x040000
  if (currentPage()) {
#else
  if (currentWidget()) {
#endif

    if ( isTabSelected()) {

#if QT_VERSION < 0x040000
      QApplication::sendPostedEvents () ;
#else
      QCoreApplication::sendPostedEvents () ;
#endif

#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4QTabWidget::paintEvent OK\n");
#endif
#if QT_VERSION < 0x040000
      QString text = label (currentPageIndex()); 
#else
      QString text = tabText (currentIndex()); 
#endif

      if (lastCreated == -1) {
        QString paramSelect = QString("/vis/viewer/select ")+text;
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if(UI != NULL)  {
#if QT_VERSION < 0x040000
          UI->ApplyCommand(paramSelect.ascii());
#else
          UI->ApplyCommand(paramSelect.toStdString().c_str());
#endif
        }
      } else {
        lastCreated = -1;
      }
      setTabSelected(false);
    }
  }
}

#endif
