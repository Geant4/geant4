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
// $Id: G4UIQt.cc,v 1.30 2010-01-06 14:13:08 lgarnier Exp $
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
#include <qresizeevent>
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
:fG4cout("")
,fHelpTBWidget(NULL)
,fHistoryTBWidget(NULL)
,fCoutTBWidget(NULL)
,fVisParametersTBWidget(NULL)
,fViewComponentsTBWidget(NULL)
,fTabWidget(NULL)
,fCoutText("Output")
{

  G4Qt* interactorManager = G4Qt::getInstance (argc,argv,(char*)"Qt");
  if (!(QApplication*)interactorManager->GetMainInteractor()) {
    G4cout        << "G4UIQt : Unable to init Qt. Aborted" << G4endl;
  }
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

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
  printf("G4UIQt::Initialise after main window creation UImanager:%d G4UIQt(this) : %d+++++++++++\n",UI,this);
#endif

  QWidget *mainWidget = new QWidget(fMainWindow);
#if QT_VERSION < 0x040000
  QSplitter *myVSplitter = new QSplitter(Qt::Horizontal,fMainWindow);
  fToolBox = new QToolBox(Qt::Horizontal,myVSplitter);
#else
  QSplitter *myVSplitter = new QSplitter(Qt::Horizontal,fMainWindow);
  fToolBox = new QToolBox(mainWidget);
#endif

  // Set layouts

  // Add a empty tabwidget
  fTabWidget = new QTabWidget(myVSplitter);

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
#if QT_VERSION >= 0x040000
  QVBoxLayout *mainLayout = new QVBoxLayout();
#else
  QVBoxLayout *mainLayout = new QVBoxLayout(mainWidget);
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
  fToolBox->addItem(fVisParametersTBWidget,"Vis parameters");
  fToolBox->addItem(fViewComponentsTBWidget,"Viewer components");
  fToolBox->addItem(fHelpTBWidget,"Help");
  fToolBox->addItem(fCoutTBWidget,"Cout");
  fToolBox->addItem(fHistoryTBWidget,"History");



  fToolBox->setSizePolicy (QSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed));
  fTabWidget->setSizePolicy (QSizePolicy(QSizePolicy::Maximum,QSizePolicy::Maximum));

  QSizePolicy policy = fTabWidget->sizePolicy();
  policy.setHorizontalStretch(1);
  policy.setVerticalStretch(1);
  fTabWidget->setSizePolicy(policy);

  myVSplitter->addWidget(fToolBox);
  myVSplitter->addWidget(fTabWidget);

#if QT_VERSION >= 0x040000
  commandLineWidget->setLayout(layoutCommandLine);
#endif
  commandLineWidget->setSizePolicy (QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));
  mainLayout->addWidget(myVSplitter,1);
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

  connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));

  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP

  interactorManager->SetG4UI(this);

#if QT_VERSION < 0x040000
  fMainWindow->setCaption( tr( "G4UI Session" ));
  fMainWindow->resize(900,600); 
  fMainWindow->move(50,100);
#else
  fMainWindow->setWindowTitle( tr("G4UI Session") ); 
  fMainWindow->resize(900,600); 
  fMainWindow->move(QPoint(50,100));
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
  QVBoxLayout *vLayout = new QVBoxLayout(fHelpTBWidget);
  QSplitter *splitter = new QSplitter(Qt::Horizontal,fHelpTBWidget);
#else
  QWidget *helpWidget = new QWidget();
  QHBoxLayout *helpLayout = new QHBoxLayout();
  QVBoxLayout *vLayout = new QVBoxLayout();
  QSplitter *splitter = new QSplitter(Qt::Horizontal);
#endif
  helpLine = new QLineEdit(fHelpTBWidget);
  helpLayout->addWidget(new QLabel("Search :",helpWidget));
  helpLayout->addWidget(helpLine);
#if QT_VERSION < 0x040000
  connect( helpLine, SIGNAL( returnPressed () ), this, SLOT( LookForHelpStringCallback() ) );
#else
  connect( helpLine, SIGNAL( editingFinished () ), this, SLOT( LookForHelpStringCallback() ) );
#endif
  
  // the help tree
#if QT_VERSION < 0x040000
  fHelpTreeWidget = new QListView(splitter);
#else
  fHelpTreeWidget = new QTreeWidget();
#endif
  fHelpTreeWidget = CreateHelpTree();
  
  fHelpArea = new QTextEdit(splitter);
  fHelpArea->setReadOnly(true);
  
  // Set layouts
  
#if QT_VERSION >= 0x040000
  if (fHelpTreeWidget)
    splitter->addWidget(fHelpTreeWidget);
  splitter->addWidget(fHelpArea);
#endif
  
  
#if QT_VERSION >= 0x040000
  vLayout->addWidget(helpWidget);
  vLayout->addWidget(splitter,1);
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
  splitter->setSizes(list);
  
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
void G4UIQt::AddTabWidget(
 QWidget* aWidget
,QString name
,int sizeX
,int sizeY
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::AddTabWidget %d %d\n",sizeX, sizeY);
#endif
  if (!aWidget) return;
  
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::AddTabWidget ADD %d %d + %d %d---------------------------------------------------\n",sizeX, sizeY,sizeX-fTabWidget->width(),sizeY-fTabWidget->height());
#endif
  
#if QT_VERSION < 0x040000
  aWidget->reparent(fTabWidget,0,QPoint(0,0));  
#else
  aWidget->setParent(fTabWidget);
#endif
    
  fTabWidget->insertTab(fTabWidget->count()-1,aWidget,name);
  
  fMainWindow->resize(fMainWindow->width()+sizeX-fTabWidget->width(),
                      fMainWindow->height()+sizeY-fTabWidget->height());
}


void G4UIQt::UpdateTabWidget(int tabNumber) {
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::UpdateTabWidget\n");
#endif
  if (  fTabWidget == NULL) {
    fTabWidget = new QTabWidget;
  }
  
  QString text = fTabWidget->tabText (tabNumber); 
  if ( text != fCoutText) {

    // Then we have to do /vis/viewer/select on the selelected tab
    // And resize selected tab
    fTabWidget->currentWidget()->resize(fTabWidget->currentWidget()->width(),fTabWidget->currentWidget()->height());

    QString paramSelect = QString("/vis/viewer/select ")+text;
    QString paramFlush = QString("/vis/viewer/flush ")+text;
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if(UI!=NULL)  {
      UI->ApplyCommand(paramSelect.toStdString().c_str());
      UI->ApplyCommand(paramFlush.toStdString().c_str());
    }
    fTabWidget->currentWidget()->repaint();
  }
  if (fTabWidget->count() == 1) {
    connect(fTabWidget, SIGNAL(currentChanged ( int ) ), SLOT(UpdateTabWidget(int)));
    //  connect(fTabWidget, SIGNAL(resizeEvent (  QResizeEvent* ) ), SLOT(ResizeTabWidget( QResizeEvent*)));
  }
}


/** Send resize event to all tabs
 */
void G4UIQt::ResizeTabWidget( QResizeEvent* e) {
  for (int a=0;a<fTabWidget->count() ;a++) {
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

  G4Qt* interactorManager = G4Qt::getInstance ();

#if QT_VERSION >= 0x040000
#if QT_VERSION >= 0x040200
  fMainWindow->setVisible(true);
#else
  fMainWindow->show();
#endif
#else
  fMainWindow->show();
#endif
  Prompt("Session :");
  exitSession = false;


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
 
  QStringList result = newStr.filter(fCoutFilter->text());

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
 
 QStringList result = newStr.filter(fCoutFilter->text());

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
  QAction *action = new QAction(QString(aLabel),QString(aLabel),QKeySequence::QKeySequence (),signalMapper, SLOT(map()));
  action->addTo(parent);
 connect(action,SIGNAL(activated()),signalMapper,SLOT(map()));

#elif QT_VERSION < 0x040000
  QAction *action = new QAction(QString(aLabel),QKeySequence::QKeySequence (),signalMapper, SLOT(map()));
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
  printf("G4UIQt::ActivateCommand found : %d \n",targetCom.data());
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

#if QT_VERSION < 0x040000
QListView * G4UIQt::CreateHelpTree()
#else
QTreeWidget * G4UIQt::CreateHelpTree()
#endif
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return NULL;
  G4UIcommandTree * treeTop = UI->GetTree();


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

  G4int treeSize = treeTop->GetTreeEntry();
#if QT_VERSION < 0x040000
  QListViewItem * newItem;
#else
  QTreeWidgetItem * newItem;
#endif
  for (int a=0;a<treeSize;a++) {
    // Creating new item

#if QT_VERSION < 0x040000
    newItem = new QListViewItem(fHelpTreeWidget);
    newItem->setText(0,QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).simplifyWhiteSpace());
#else
    newItem = new QTreeWidgetItem(fHelpTreeWidget);
    newItem->setText(0,QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).trimmed());
#endif


    // look for childs
    CreateChildTree(newItem,treeTop->GetTree(a+1));
    //      items.append(newItem);
  }


#if QT_VERSION < 0x040000
  connect(fHelpTreeWidget, SIGNAL(selectionChanged ()),this, SLOT(HelpTreeClicCallback()));  
  connect(fHelpTreeWidget, SIGNAL(doubleClicked (QListViewItem*)),this, SLOT(HelpTreeDoubleClicCallback()));
#else
  connect(fHelpTreeWidget, SIGNAL(itemSelectionChanged ()),this, SLOT(HelpTreeClicCallback()));  
  connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked (QTreeWidgetItem*,int)),this, SLOT(HelpTreeDoubleClicCallback()));  
#endif

  return fHelpTreeWidget;
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

  // Get the Sub directories
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {

#if QT_VERSION < 0x040000
      newItem = new QListViewItem(aParent);
      newItem->setText(0,QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).simplifyWhiteSpace());

#else
    newItem = new QTreeWidgetItem(aParent);
    newItem->setText(0,QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).trimmed());
#endif

    CreateChildTree(newItem,aCommandTree->GetTree(a+1));
  }



  // Get the Commands

  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    QStringList stringList;
#if QT_VERSION < 0x040000
    newItem = new QListViewItem(aParent);
    newItem->setText(0,QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).simplifyWhiteSpace());
    newItem->setOpen(false);

#else
    newItem = new QTreeWidgetItem(aParent);
    newItem->setText(0,QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed());
#if QT_VERSION < 0x040202
   	fHelpTreeWidget->setItemExpanded(newItem,false); 
#else
    newItem->setExpanded(false);
#endif
#endif
      
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

  if (aParent->text(0) == aCommand)
    return aParent;
  
#if QT_VERSION < 0x040000
  QListViewItem * tmp = NULL;
  QListViewItem* tmpItem = aParent->firstChild();
    while (tmpItem != 0) {
      if (!tmp)
        tmp = FindTreeItem(tmpItem,aCommand);
      tmpItem = tmpItem->nextSibling();
    }

#else
  QTreeWidgetItem * tmp = NULL;
  for (int a=0;a<aParent->childCount();a++) {
    if (!tmp)
      tmp = FindTreeItem(aParent->child(a),aCommand);
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



  std::string itemText;
#if QT_VERSION < 0x040000
  itemText = std::string(item->text(0).ascii());
#else
  itemText = std::string(item->text(0).toStdString());
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
  fCommandArea->setText(item->text(0));
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
const QString & text) {

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::CoutFilterCallback : %s\n",text.toStdString().c_str());
#endif
 
  QStringList result = fG4cout.filter(text);

  fCoutTBTextArea->setText(result.join("\n"));
  fCoutTBTextArea->repaint();
#if QT_VERSION < 0x040000
  fCoutTBTextArea->verticalScrollBar()->setValue(fCoutTBTextArea->verticalScrollBar()->maxValue());
#else
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
#endif
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIQt::CoutFilterCallback textSize:%d resSize:%d \n",fG4cout.join("").size(),result.join("").size());
#endif

 }

/**   Callback called when user give a new string to look for<br>
   Display a list of matching commands descriptions. If no string is set,
   will display the complete help tree
*/
void G4UIQt::LookForHelpStringCallback(
)
{
  QString searchText = helpLine->text();

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
    CreateHelpTree();
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
    newItem->setText(0,i.data().simplifyWhiteSpace());
#else
    newItem = new QTreeWidgetItem(fHelpTreeWidget);
    newItem->setText(0,i.value().trimmed());
#endif
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
#endif
