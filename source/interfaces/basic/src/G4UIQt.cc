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
// $Id: G4UIQt.cc,v 1.53 2010-11-02 15:38:51 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include <qmenu.h>
#include <qlistwidget.h>
#include <qtreewidget.h>



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
  Q_FOREACH (QWidget *widget, QApplication::allWidgets()) {
    if ((found== false) && (widget->inherits("QMainWindow"))) {
      found = true;
    }
  }

  if (found) {
    G4cout        << "G4UIQt : Found an external App with a QMainWindow already defined. Aborted" << G4endl;
    return ;
  }
  fMainWindow = new QMainWindow();

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::Initialise after main window creation +++++++++++\n");
#endif

  QWidget *mainWidget = new QWidget(fMainWindow);
  fMyVSplitter = new QSplitter(Qt::Horizontal,fMainWindow);
  fToolBox = new QToolBox();

  // Set layouts

  QWidget* commandLineWidget = new QWidget(mainWidget);
  QVBoxLayout *layoutCommandLine = new QVBoxLayout();

  // fill them

  fCommandLabel = new QLabel("",commandLineWidget);

  fCommandArea = new QLineEdit(commandLineWidget);
  fCommandArea->installEventFilter(this);
  fCommandArea->activateWindow();

  fCommandArea->setFocusPolicy ( Qt::StrongFocus );
  fCommandArea->setFocus(Qt::TabFocusReason);



  layoutCommandLine->addWidget(fCommandLabel);
  layoutCommandLine->addWidget(fCommandArea);
  QVBoxLayout *mainLayout;
  mainLayout = new QVBoxLayout();

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

  fEmptyViewerTabLabel = new QLabel("         If you want to have a Viewer, please use /vis/open commands. ");

  // Only at creation. Will be set visible when sessionStart();
 #if QT_VERSION < 0x040200
  fEmptyViewerTabLabel->hide();
 #else
  fEmptyViewerTabLabel->setVisible(false);
 #endif


  fMyVSplitter->addWidget(fToolBox);
  fMyVSplitter->addWidget(fEmptyViewerTabLabel);

  commandLineWidget->setLayout(layoutCommandLine);
  commandLineWidget->setSizePolicy (QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));
  mainLayout->addWidget(fMyVSplitter,1);
  mainLayout->addWidget(commandLineWidget);

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt :: 5\n");
#endif

  mainWidget->setLayout(mainLayout);

  fMainWindow->setCentralWidget(mainWidget);


  // Add a quit subMenu
  QMenu *fileMenu = fMainWindow->menuBar()->addMenu("File");
  fileMenu->addAction("Quit", this, SLOT(ExitSession()));


  AddInteractor ("file",(G4Interactor)fileMenu);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::G4UIQt :: 6\n");
#endif

  // Connect signal
  connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));
  connect(fToolBox, SIGNAL(currentChanged(int)), SLOT(ToolBoxActivated(int)));

  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP

  fMainWindow->setWindowTitle( tr("G4UI Session") ); 
  fMainWindow->resize(900,600); 
  fMainWindow->move(QPoint(50,100));

  // Set not visible until session start
 #if QT_VERSION < 0x040200
  fMainWindow->hide();
 #else
  fMainWindow->setVisible(false);
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

  QVBoxLayout *layoutHistoryTB = new QVBoxLayout();
  fHistoryTBTableList = new QListWidget();
  fHistoryTBTableList->setSelectionMode(QAbstractItemView::SingleSelection);
  connect(fHistoryTBTableList, SIGNAL(itemSelectionChanged()), SLOT(CommandHistoryCallback()));
  fHistoryTBTableList->installEventFilter(this);

  layoutHistoryTB->addWidget(fHistoryTBTableList);

  fHistoryTBWidget->setLayout(layoutHistoryTB);
}

/** Create the Help ToolBox Widget
 */
void G4UIQt::CreateHelpTBWidget(
) 
{

  
  QWidget *helpWidget = new QWidget();
  QHBoxLayout *helpLayout = new QHBoxLayout();
  QVBoxLayout *vLayout = new QVBoxLayout();
  fHelpVSplitter = new QSplitter(Qt::Horizontal);
  fHelpLine = new QLineEdit(fHelpTBWidget);
  helpLayout->addWidget(new QLabel("Search :",helpWidget));
  helpLayout->addWidget(fHelpLine);
  connect( fHelpLine, SIGNAL( editingFinished () ), this, SLOT( LookForHelpStringCallback() ) );
  
  // Create Help tree
  FillHelpTree();
  
  fHelpArea = new QTextEdit(fHelpVSplitter);
  fHelpArea->setReadOnly(true);
  
  // Set layouts
  
  if (fHelpTreeWidget) {
    fHelpVSplitter->addWidget(fHelpTreeWidget);
  }
  fHelpVSplitter->addWidget(fHelpArea);
  
  
  vLayout->addWidget(helpWidget);
  vLayout->addWidget(fHelpVSplitter,1);
  
  fHelpTBWidget->setMinimumSize(50,50);
  fHelpTBWidget->setSizePolicy (QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));
  // set the splitter size
  QList<int> list;
  list.append( 50 );
  list.append( 50 );
  fHelpVSplitter->setSizes(list);
  
  helpWidget->setLayout(helpLayout);
  fHelpTBWidget->setLayout(vLayout);
}


/** Create the Cout ToolBox Widget
 */
void G4UIQt::CreateCoutTBWidget(
) 
{
  QVBoxLayout *layoutCoutTB = new QVBoxLayout();

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

  fCoutTBWidget->setLayout(layoutCoutTB);
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
#if QT_VERSION < 0x040500
#else
    fTabWidget->setTabsClosable (true); 
#endif
    
#if QT_VERSION < 0x040200
#else
    fTabWidget->setUsesScrollButtons (true);
#endif
    
    fTabWidget->setSizePolicy (QSizePolicy(QSizePolicy::Maximum,QSizePolicy::Maximum));
    
    QSizePolicy policy = fTabWidget->sizePolicy();
    policy.setHorizontalStretch(1);
    policy.setVerticalStretch(1);
    fTabWidget->setSizePolicy(policy);
    
#if QT_VERSION < 0x040500
#else
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
  if (fEmptyViewerTabLabel != NULL) {
    if ( fMyVSplitter->indexOf(fEmptyViewerTabLabel) != -1) {
      
      fEmptyViewerTabLabel->hide();
      fEmptyViewerTabLabel->setParent(NULL);
      delete fEmptyViewerTabLabel;
      fEmptyViewerTabLabel = NULL;
      
      fMyVSplitter->addWidget(fTabWidget);
      
      aWidget->setParent(fTabWidget);
    }
  }

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::AddTabWidget ADD %d %d + %d %d---------------------------------------------------\n",sizeX, sizeY,sizeX-fTabWidget->width(),sizeY-fTabWidget->height());
#endif

  if (fMainWindow->isVisible()) {

    // get the size of the tabbar
    int tabBarX = 0;
    int tabBarY = 0;
    if (fTabWidget->count() >0) {
      tabBarX = fTabWidget->width()-fTabWidget->widget(0)->width();
      tabBarY = fTabWidget->height()-fTabWidget->widget(0)->height();
    }

    fMainWindow->resize(tabBarX+fMainWindow->width()+sizeX-fTabWidget->width(),tabBarY+fMainWindow->height()+sizeY-fTabWidget->height());
  }

  // Problems with resize. The widgets are not realy drawn at this step,
  // then we have to force them on order to check the size

  fTabWidget->insertTab(fTabWidget->count(),aWidget,name);
  
  fTabWidget->setCurrentIndex(fTabWidget->count()-1);

  // Set visible
 #if QT_VERSION < 0x040200
   fTabWidget->setLastTabCreated(fTabWidget->currentIndex());
 #else
   fTabWidget->setLastTabCreated(fTabWidget->currentIndex());
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

  fTabWidget->setCurrentIndex(tabNumber);

  // Send this signal to unblock graphic updates !
  fTabWidget->setTabSelected(false);

 #if QT_VERSION < 0x040200
  fTabWidget->show();
 #else
  fTabWidget->setVisible(true);
 #endif

  // This will send a paintEvent to OGL Viewers
  fTabWidget->setTabSelected(true);

  QCoreApplication::sendPostedEvents () ;

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
    fTabWidget->widget(a)->resize(e->size());
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

  #if QT_VERSION < 0x040200
      fMainWindow->show();
  #else
      fMainWindow->setVisible(true);
  #endif
  // get the size of the tabbar
  int tabBarX = 0;
  int tabBarY = 0;

  if (fTabWidget != NULL) {
    tabBarX = -fTabWidget->widget(0)->width();
    tabBarY = -fTabWidget->widget(0)->height();
  }
  fMainWindow->resize(tabBarX+fMainWindow->width()+fLastQTabSizeX,tabBarY+fMainWindow->height()+fLastQTabSizeY);

  QCoreApplication::sendPostedEvents () ;

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
  newStr = QStringList(QString((char*)aString.data()).trimmed());
  fG4cout += newStr;
 
  QStringList result = newStr.filter(fCoutFilter->text());

  if (result.join("\n").isEmpty()) {
    return 0;
  }
  fCoutTBTextArea->append(result.join("\n"));
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
 G4String aString
)
{
  if (!aString) return 0;

  QStringList newStr;

  // Add to stringList
  newStr = QStringList(QString((char*)aString.data()).trimmed());
  fG4cout += newStr;
 
  QStringList result = newStr.filter(fCoutFilter->text());

  QColor previousColor = fCoutTBTextArea->textColor();
  fCoutTBTextArea->setTextColor(Qt::red);
  fCoutTBTextArea->append(result.join("\n"));
  fCoutTBTextArea->setTextColor(previousColor);
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
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

  QMenu *fileMenu = new QMenu(aLabel);
  fMainWindow->menuBar()->insertMenu(fMainWindow->menuBar()->actions().last(),fileMenu); 

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

  QMenu *parent = (QMenu*)GetInteractor(aMenu);

  if(parent==NULL) return;
  
  QSignalMapper *signalMapper = new QSignalMapper(this);
  QAction *action = parent->addAction(aLabel, signalMapper, SLOT(map()));

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

  fToolBox->setCurrentWidget(fHelpTBWidget);
}



/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIQt::InitHelpTree()
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
    InitHelpTree();
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
      
      newItem = new QTreeWidgetItem(fHelpTreeWidget);
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
void G4UIQt::CreateChildTree(
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
      newItem = new QTreeWidgetItem(aParent);
      newItem->setText(0,GetShortCommandPath(commandText));
    }
    CreateChildTree(newItem,aCommandTree->GetTree(a+1));
  }



  // Get the Commands

  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    QStringList stringList;
    commandText = QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed();

    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
      newItem = new QTreeWidgetItem(aParent);
      newItem->setText(0,GetShortCommandPath(commandText));
#if QT_VERSION < 0x040202
      fHelpTreeWidget->setItemExpanded(newItem,false); 
#else
      newItem->setExpanded(false);
#endif
    }
  }
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
  G4String command (fCommandArea->text().toStdString().c_str());
  if (fCommandArea->text().trimmed() != "") {
    fHistoryTBTableList->addItem(fCommandArea->text());
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
  G4String ss = G4String(aCommand.toStdString().c_str());
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
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIQt::CommandHistoryCallback change text\n");
#endif
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
    newItem = new QTreeWidgetItem(fHelpTreeWidget);
    QString commandStr = i.value().trimmed();

    if (commandStr.indexOf("/") == 0) {
      commandStr = commandStr.right(commandStr.size()-1);
    }
      
    newItem->setText(0,commandStr);
    newItem->setText(1,progressStr);
    
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


  
#if QT_VERSION < 0x040500
void G4UIQt::TabCloseCallback(int){
#else
void G4UIQt::TabCloseCallback(int a){
#endif
#if QT_VERSION < 0x040500
#else
  QWidget* temp = fTabWidget->widget(a);
  fTabWidget->removeTab (a);

  delete temp;

  if (fTabWidget->count() == 0) {
    if (fEmptyViewerTabLabel == NULL) {
      fEmptyViewerTabLabel = new QLabel("         If you want to have a Viewer, please use /vis/open commands. ");
    }

    fMyVSplitter->addWidget(fEmptyViewerTabLabel);
    fMyVSplitter->show();
    fEmptyViewerTabLabel->show();
    fTabWidget->setParent(0);
    fTabWidget->setVisible(false);
    delete fTabWidget;
    fTabWidget = NULL;
  }
#endif
}


void G4UIQt::ToolBoxActivated(int a){

  if (fToolBox->widget(a) == fHelpTBWidget) {
    // Rebuild the help tree
    FillHelpTree();
  }
}

void G4QTabWidget::paintEvent(
QPaintEvent *
)
{

  if (currentWidget()) {

    if ( isTabSelected()) {

      QCoreApplication::sendPostedEvents () ;

#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4QTabWidget::paintEvent OK\n");
#endif
      QString text = tabText (currentIndex()); 

      if (lastCreated == -1) {
        QString paramSelect = QString("/vis/viewer/select ")+text;
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if(UI != NULL)  {
          UI->ApplyCommand(paramSelect.toStdString().c_str());
        }
      } else {
        lastCreated = -1;
      }
      setTabSelected(false);
    }
  }
}

#endif
