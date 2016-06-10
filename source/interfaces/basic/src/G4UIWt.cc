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
// $Id: G4UIWt.cc,v 1.53 2010-11-02 15:38:51 lgarnier Exp $
//
// L. Garnier

#ifdef G4UI_BUILD_WT_SESSION


#include <string.h>

#include "G4UIWt.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"

#include "G4Wt.hh"

#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WTextArea>
#include <Wt/WMessageBox>
#include <Wt/WLayout>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WComboBox>
#include <Wt/WVBoxLayout>
#include <Wt/WHBoxLayout>
#include <Wt/WGridLayout>
#include <Wt/WPanel>
#include <Wt/WSelectionBox>
#include <Wt/WLength>

#define G4DEBUG_INTERFACES_BASIC 1

// Pourquoi Static et non  variables de classe ?
static G4bool exitSession = true;
static G4bool exitPause = true;

/**   Build a Wt window with a menubar, output area and promt area<br> 
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
G4UIWt::G4UIWt (
 int argc
,char** argv
)
:fMainWindow(NULL)
,fCommandLabel(NULL)
,fCommandArea(NULL)
,fCoutTBTextArea(NULL)
,fHelpArea(NULL)
,fUITabWidget(NULL)
,fG4cout()
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
,fHelpVSplitter(NULL)
,fToolbarApp(NULL)
,fToolbarUser(NULL)
,fStringSeparator("__$$$@%%###__")
,fLastOpenPath("")
,fMoveSelected(false)
,fRotateSelected(true)
,fPickSelected(false)
,fZoomInSelected(false)
,fZoomOutSelected(false)
,fExitSession(true)
,fExitPause(true)

{
  
  G4Wt* interactorManager = G4Wt::getInstance (argc,argv,(char*)"Wt");
  if (!(Wt::WApplication*)interactorManager->GetMainInteractor()) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout        << "G4UIWt : Unable to init Wt. Aborted" << G4endl;
    }
  }
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);
  if(UI!=NULL) UI->SetG4UIWindow(this);
  
  fMainWindow = new Wt::WContainerWidget(wApp->root());
  
  // resize to a big size in order to gave space for the viewer
  Wt::WVBoxLayout* mainWindowVLayout = new Wt::WVBoxLayout();
  fMainWindow->setLayout(mainWindowVLayout);
  

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::Initialise after main window creation +++++++++++\n");
#endif


  // the splitter
  fMainSplitterWidget = new Wt::WContainerWidget();
  Wt::WGridLayout* fMainSplitterWidgetLayout = new Wt::WGridLayout();


  fMainSplitterWidgetLayout->addWidget(CreateLeftSplitterWidget(),1,1);
  fMainSplitterWidgetLayout->addWidget(CreateRightSplitterWidget(),1,2);
  fMainSplitterWidgetLayout->setColumnResizable	(1,true,Wt::WLength(30,Wt::WLength::Percentage));
  
  // create vis tab widget
//  Wt::WTabWidget* tabWidget = new Wt::WTabWidget();
  

  fMainSplitterWidget->setLayout(fMainSplitterWidgetLayout);
  mainWindowVLayout->addWidget(fMainSplitterWidget);



#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt :: 5\n");
#endif


  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP

  wApp->setTitle("Wt application");
  
  // Set not visible until session start
  fMainWindow->hide();
}



G4UIWt::~G4UIWt(
) 
{ 
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::~G4UIWt Delete\n");
#endif
  G4UImanager* UI = G4UImanager::GetUIpointer();  // TO KEEP
  if(UI!=NULL) {  // TO KEEP
    UI->SetSession(NULL);  // TO KEEP
    UI->SetG4UIWindow(NULL);
    UI->SetCoutDestination(NULL);  // TO KEEP
  }
  
  if (fMainWindow!=NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::~G4UIWt DELETE fMainWindow\n");
#endif
    delete fMainWindow;
  }
}

/** Create the History ToolBox Widget
 */
Wt::WWidget* G4UIWt::CreateHistoryTBWidget(
)
{
  fHistoryTBWidget = new Wt::WPanel();
  fHistoryTBWidget->setCollapsed(true);
  fHistoryTBWidget->setCollapsible(true);
  
  fHistoryTBTableList = new Wt::WSelectionBox();
  fHistoryTBTableList->setSelectionMode(Wt::SingleSelection);
  fHistoryTBTableList->changed().connect(this,&G4UIWt::CommandHistoryCallback);
  
  fHistoryTBWidget->setCentralWidget(fHistoryTBTableList);
  return fHistoryTBWidget;
}


/** Create the Help ToolBox Widget
 */
Wt::WWidget* G4UIWt::CreateHelpTBWidget(
) 
{
  fHelpTBWidget = new Wt::WPanel();
  fHelpTBWidget->setCollapsible(true);
  
  Wt::WContainerWidget *helpWidget = new Wt::WContainerWidget();
  Wt::WHBoxLayout *helpLayout = new Wt::WHBoxLayout();
  Wt::WVBoxLayout *vLayout = new Wt::WVBoxLayout();
  
  fHelpVSplitter = new Wt::WContainerWidget();
  fHelpVSplitter->setLayout(new Wt::WHBoxLayout());
  
  Wt::WVBoxLayout* VHelpSplitterVLayout = new Wt::WVBoxLayout();
  fHelpVSplitter->setLayout(VHelpSplitterVLayout);
  
  
  fHelpLine = new Wt::WLineEdit();
  fHelpTBWidget->setCentralWidget(helpWidget);
  helpWidget->setLayout(helpLayout);
  helpWidget->layout()->addWidget(new Wt::WLabel(Wt::WString("Search :")));
  helpWidget->layout()->addWidget(fHelpLine);
  printf("*** G4UIWt::CreateHelpTBWidget, missing EnterPress connection\n");
  //  fHelpLine->enterPressed().connect(this, &G4UIWt::LookForHelpStringCallback );
  
  // Create Help tree
  FillHelpTree();
  
  fHelpArea = new Wt::WTextArea(fHelpVSplitter);
  fHelpArea->setReadOnly(true);
  
  // Set layouts
  
  if (fHelpTreeWidget) {
    fHelpVSplitter->addWidget(fHelpTreeWidget);
  }
  fHelpVSplitter->addWidget(fHelpArea);
  
  vLayout->addWidget(helpWidget);
  vLayout->addWidget(fHelpVSplitter,1);
  
  //  ((Wt::WPanel*) fHelpTBWidget)->setMinimumSize(50,50);
  //  fHelpTBWidget->setSizePolicy (WSizePolicy(WSizePolicy::Minimum,WSizePolicy::Minimum));
  // set the splitter size
  //  Wt::WList<int> list;
  //   list.append( 50 );
  //   list.append( 50 );
  //  fHelpVSplitter->setSizes(list);
  
  return fHelpTBWidget;
}


/** Create the Cout ToolBox Widget
 */
Wt::WWidget* G4UIWt::CreateCoutTBWidget(
) 
{
  fCoutTBWidget = new Wt::WPanel();
  fCoutTBWidget->setCollapsible(true);
  fCoutTBWidget->setTitle("Output");
  fCoutTBWidget->setCollapsed(false);
  
  Wt::WContainerWidget* myContainer = new Wt::WContainerWidget();
  Wt::WVBoxLayout *myContainerVLayout = new Wt::WVBoxLayout();
  myContainer->setLayout(myContainerVLayout);
  
  // Could have been created if any err ou cout *before* sessionStart()
  if (!fCoutTBTextArea) {
    fCoutTBTextArea = new Wt::WTextArea(myContainer);
  }
  fCoutFilter = new Wt::WLineEdit();
  Wt::WLabel* coutFilterLabel = new Wt::WLabel("Filter : ");
  myContainerVLayout->addWidget(fCoutTBTextArea);
  myContainerVLayout->addWidget(fCoutFilter);
  myContainerVLayout->addWidget(coutFilterLabel);
  
  Wt::WPushButton *coutTBClearButton = new Wt::WPushButton ("clear",myContainer);
  
  //  fCoutFilter->changed().connect(this, (&G4UIWt::CoutFilterCallback));
  
  // SIGNAL/SLOT connection
  // To be implemented
  
  fCoutTBTextArea->setReadOnly(true);
  
  Wt::WContainerWidget* coutButtonWidget = new Wt::WContainerWidget(myContainer);
  Wt::WHBoxLayout* layoutCoutTBButtons = new Wt::WHBoxLayout(coutButtonWidget);
  layoutCoutTBButtons->addWidget(coutTBClearButton);
  layoutCoutTBButtons->addWidget(coutFilterLabel);
  layoutCoutTBButtons->addWidget(fCoutFilter);
  
  myContainerVLayout->addWidget(coutButtonWidget);
  
  fCoutTBWidget->setCentralWidget(myContainer);

  return fCoutTBWidget;
}


/** Create the VisParameters ToolBox Widget
 */
Wt::WContainerWidget* G4UIWt::CreateVisParametersTBWidget(
) 
{
  return NULL;
}


/** Create the VisParameters ToolBox Widget
 */
Wt::WWidget* G4UIWt::CreateUITabWidget(
) 
{
  fUITabWidget = new Wt::WTabWidget();

  // the right splitter
  fUITabWidget->addTab(CreateSceneTreeComponentsTBWidget(),"Scene tree");
  fUITabWidget->addTab(CreateHelpTBWidget(),"Help");
  fUITabWidget->addTab(CreateHistoryTBWidget(),"History");
  //  fUITabWidget->setCurrentWidget(fSceneTreeComponentsTBWidget);

  // SIGNAL/SLOT connection
  // To be implemented
  
  
  return fUITabWidget;
}


Wt::WWidget* G4UIWt::CreateSceneTreeComponentsTBWidget(){

  fSceneTreeComponentsTBWidget = new Wt::WTabWidget();

  fSceneTreeComponentsTBWidget->hide();
  
  return fSceneTreeComponentsTBWidget;
}


Wt::WContainerWidget* G4UIWt::CreateLeftSplitterWidget(){
  
  fLeftSplitterWidget = new Wt::WContainerWidget();
  Wt::WVBoxLayout * layoutLeftSplitterWidget = new Wt::WVBoxLayout();
  fLeftSplitterWidget->setLayout(layoutLeftSplitterWidget);

  layoutLeftSplitterWidget->addWidget(CreateUITabWidget());

  return fLeftSplitterWidget;
}


Wt::WContainerWidget* G4UIWt::CreateRightSplitterWidget(){

  fRightSplitterWidget = new Wt::WContainerWidget();

  // Set layouts
  Wt::WVBoxLayout* VSplitterVLayout = new Wt::WVBoxLayout();
  fRightSplitterWidget->setLayout(VSplitterVLayout);
  Wt::WContainerWidget* commandLineWidget = new Wt::WContainerWidget();
  
  Wt::WHBoxLayout *layoutCommandLine = new Wt::WHBoxLayout();
  commandLineWidget->setLayout(layoutCommandLine);

  // fill them

  fCommandLabel = new Wt::WLabel("",commandLineWidget);
  fCommandArea = new Wt::WLineEdit(commandLineWidget);
  fCommandArea->setToolTip("Apply command");


  layoutCommandLine->addWidget(fCommandLabel);
  layoutCommandLine->addWidget(fCommandArea);

  
  fViewerTabWidget = new G4WTabWidget(fRightSplitterWidget);
  
  fViewerTabWidget->tabClosed().connect(this, &G4UIWt::TabCloseCallback);
  
  fViewerTabWidget->currentChanged().connect(this,&G4UIWt::CurrentChangedTabWidgetCallback);

  
  fEmptyViewerTabLabel = new Wt::WLabel("         If you want to have a Viewer, please use /vis/open commands. ");

  // UI Specific fonctions
  fRightSplitterWidget->layout()->addWidget(fViewerTabWidget);
  fRightSplitterWidget->layout()->addWidget(fEmptyViewerTabLabel);
  fRightSplitterWidget->layout()->addWidget(CreateCoutTBWidget());
  fRightSplitterWidget->layout()->addWidget(commandLineWidget);

  commandLineWidget->setMinimumSize(50,50);

  // Connect signal
  fCommandArea->enterPressed().connect(this,&G4UIWt::CommandEnteredCallback);

  return fRightSplitterWidget;
}


/** Get the ViewerComponents ToolBox Widget
 */
Wt::WTabWidget* G4UIWt::GetSceneTreeComponentsTBWidget(
)
{
  return fSceneTreeComponentsTBWidget;
}


/**   Add a new tab widget.
 Create the tab if it was not done
 */
bool G4UIWt::AddTabWidget(
                          Wt::WWidget* aWidget
                          ,Wt::WString name
                          ,int /* width */
                          ,int /* height */
                          )
{
/*
 if (fViewerTabWidget == NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget +++++\n");
#endif
 
    fViewerTabWidget = new G4WTabWidget(fRightSplitterWidget);
    
    fViewerTabWidget->tabClosed().connect(this, &G4UIWt::TabCloseCallback);
    
    fViewerTabWidget->currentChanged().connect(this,&G4UIWt::CurrentChangedTabWidgetCallback);
    fRightSplitterWidget->layout()->addWidget(fViewerTabWidget);
//    fRightSplitterWidget->layout()->addWidget(fViewerTabWidget);
    
  }
*/
  if (!fViewerTabWidget->isVisible() ) {
    if ( fRightSplitterWidget->isVisible()) {
      fRightSplitterWidget->setHidden(true);
      fEmptyViewerTabLabel->setHidden(true);
    }
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget +++++\n");
#endif
    fViewerTabWidget->show();
    fViewerTabWidget = new G4WTabWidget(fRightSplitterWidget);
    
    fViewerTabWidget->tabClosed().connect(this, &G4UIWt::TabCloseCallback);
    
    fViewerTabWidget->currentChanged().connect(this,&G4UIWt::CurrentChangedTabWidgetCallback);
    fRightSplitterWidget->layout()->addWidget(fViewerTabWidget);
    //    fRightSplitterWidget->layout()->addWidget(fViewerTabWidget);
    
  }


  if (!aWidget) {
    return false;
  }

  // Remove Wt::WLabel

  // L.Garnier 26/05/2010 : not exactly the same in qt3. Could cause some
  // troubles
  

  // Problems with resize. The widgets are not realy drawn at this step,
  // then we have to force them on order to check the size

//  aWidget->resize(width,height);
//  aWidget->resize(600,600);
//  fViewerTabWidget->resize(620,621);
  fViewerTabWidget->addTab(new Wt::WLabel("test"),"test1");
  fViewerTabWidget->addTab(aWidget,name);
  fViewerTabWidget->addTab(new Wt::WLabel("test"),"test2");

  // Change Color
  aWidget->decorationStyle().setBackgroundColor (Wt::WColor(245,245,245));

  fViewerTabWidget->setCurrentIndex(fViewerTabWidget->count()-1);

  // Set visible
  fViewerTabWidget->setLastTabCreated(fViewerTabWidget->currentIndex());

  return true;
}



void G4UIWt::CurrentChangedTabWidgetCallback(int tabNumber) {
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CurrentChangedTabWidget %d\n",tabNumber);
#endif
  if ( fViewerTabWidget == NULL) {
    fViewerTabWidget = new G4WTabWidget;
  }
  
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CurrentChangedTabWidget CALL REPAINT tabGL\n");
#endif

  fViewerTabWidget->setCurrentIndex(tabNumber);

  // Send this signal to unblock graphic updates !
  fViewerTabWidget->setTabSelected(false);

  fViewerTabWidget->show();

  // This will send a paintEvent to OGL Viewers
  fViewerTabWidget->setTabSelected(true);


  Wt::WString text = fViewerTabWidget->tabText (tabNumber);

  Wt::WString paramSelect = Wt::WString("/vis/viewer/select ")+text;
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI != NULL)  {
    UI->ApplyCommand(paramSelect.toUTF8().c_str());
  }
}


/**   Start the Wt main loop
 */
G4UIsession* G4UIWt::SessionStart (
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 1\n");
#endif
  
  G4Wt* interactorManager = G4Wt::getInstance ();
  
  Prompt("Session :");
  exitSession = false;
  
  fMainWindow->show();

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart2\n");
#endif
  interactorManager->DisableSecondaryLoop (); // TO KEEP
//  if ((Wt::WApplication*)interactorManager->GetMainInteractor())
//    ((Wt::WApplication*)interactorManager->GetMainInteractor())->exec();
  
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
void G4UIWt::Prompt (
                     G4String aPrompt
                     )
{
  if (!aPrompt) return;
  if (!fCommandLabel) return;
  
  fCommandLabel->setText((char*)aPrompt.data());
}



void G4UIWt::SessionTerminate (
)
{
  G4Wt* interactorManager = G4Wt::getInstance ();
  
  if ((Wt::WApplication*)interactorManager->GetMainInteractor()) {
    // What to do here ?
  }
}



/**
 Begin the secondary loop
 @param a_prompt : label to display as the prompt label
 */
void G4UIWt::SecondaryLoop (
                            G4String aPrompt
                            )
{
  if (!aPrompt) return;
  
  G4Wt* interactorManager = G4Wt::getInstance (); // TO KEEP ?
  Prompt(aPrompt); // TO KEEP
  fExitPause = false; // TO KEEP
  void* eventTmp; // TO KEEP
  while((eventTmp = interactorManager->GetEvent())!=NULL) {  // TO KEEP
    interactorManager->DispatchEvent(eventTmp); // TO KEEP
    if(fExitPause==true) break; // TO KEEP
  } // TO KEEP
  Prompt("Session :"); // TO KEEP
}


/**   Callback call when "click on a menu entry.<br>
 Send the associated command to geant4
 */
void G4UIWt::CommandEnteredCallback (
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CommandEnteredCallback\n");
#endif
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CommandEnteredCallback 1\n");
#endif
  G4String command (fCommandArea->text().toUTF8());
 if (fCommandArea->text() != "") {

   printf("*** G4UIWt::CommandEnteredCallback, missing update on history \n");
   /*  FIXME
    fHistoryTBTableList->addItem(fCommandArea->text());
    fHistoryTBTableList->clearSelection();
    fHistoryTBTableList->setCurrentIndex(0);
*/
   fCommandArea->setText("");
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandEnteredCallback 2\n");
#endif
    
    G4Wt* interactorManager = G4Wt::getInstance ();
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandEnteredCallback 3\n");
#endif
    if (interactorManager) {
      interactorManager->FlushAndWaitExecution();
    }
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandEnteredCallback 4\n");
#endif
    if (command(0,4) != "help") {
      ApplyShellCommand (command,exitSession,exitPause);
    } else {
      ActivateCommand(command);
    }
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandEnteredCallback 5\n");
#endif
    // Rebuild help tree
    FillHelpTree();
    
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandEnteredCallback 6\n");
#endif
    if(exitSession==true)
      SessionTerminate();
  }
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CommandEnteredCallback END\n");
#endif

 }


/**
 Add a new button to a menu
 @param aMenu : parent menu
 @param aLabel : label to display
 @param aCommand : command to execute as a callback
 */
void G4UIWt::AddButton (
                        const char* aMenu
                        ,const char* aLabel
                        ,const char* aCommand
                        )
{
  if(aMenu==NULL) return; // TO KEEP
  if(aLabel==NULL) return; // TO KEEP
  if(aCommand==NULL) return; // TO KEEP
  
  Wt::WMenu *parentTmp = (Wt::WMenu*)GetInteractor(aMenu);
  
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
  
  if(treeTop->FindPath(aCommand) == NULL) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4int verbose = UImanager->GetVerboseLevel();
    
    if (verbose >= 2) {
      G4cout << "Warning: command '"<< aCommand <<"' does not exist, please define it before using it."<< G4endl;
    }
  }
  
  printf("*** G4UIWt::AddButton, missing connection on menu callback \n");
/* FIXME
 QSignalMapper *signalMapper = new QSignalMapper(this);
  QAction *action = parentTmp->addAction(aLabel, signalMapper, SLOT(map()));
  
  connect(signalMapper, SIGNAL(mapped(const char*)),this, SLOT(ButtonCallback(const char*)));
  signalMapper->setMapping(action, aCommand);
*/
}


/**
 special case for the "open" icon. It will open a file selector and map the return file to the given command.
 */
void G4UIWt::AddIcon(const char* aLabel, const char* aIconFile, const char* aCommand, const char* /* aFileName */ ){
  if(aLabel==NULL) return; // TO KEEP
  // special case, aCommand could be NULL if aIconFile is not user_icon
  if (aCommand==NULL) {
    if (std::string(aIconFile) == "user_icon") {
      return; // TO KEEP
    }
  }
  printf("*** G4UIWt::AddIcon, missing icon creation \n");
/*
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
      "  @*{****&@@%$    $@-*-&*****+  ",
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
  Wt::WToolBar *currentToolbar = NULL;
  if (userToolBar) {
    if (fToolbarUser == NULL) {
      fToolbarUser = new Wt::WToolBar();
      fToolbarUser->setIconSize (QSize(20,20));
      fMainWindow->addToolBar(Wt::TopToolBarArea, fToolbarUser);
    }
    currentToolbar = fToolbarUser;
  } else {
    if (fToolbarApp == NULL) {
      fToolbarApp = new Wt::WToolBar();
      fToolbarApp->setIconSize (QSize(20,20));
      fMainWindow->addToolBar(Wt::TopToolBarArea, fToolbarApp);
    }
    currentToolbar = fToolbarApp;
  }
  
  QSignalMapper *signalMapper = new QSignalMapper(this);
  QAction *action = currentToolbar->addAction(pix,aLabel, signalMapper, SLOT(map()));
  
  
  // special cases :"open"
  if (std::string(aIconFile) == "open") {
    connect(signalMapper, SIGNAL(mapped(const Wt::WString &)),this, SLOT(OpenIconCallback(const Wt::WString &)));
    Wt::WString txt = aCommand + fStringSeparator + aLabel;
    signalMapper->setMapping(action, Wt::WString(txt));
    
    // special cases :"save"
  } else if (std::string(aIconFile) == "save") {
    connect(signalMapper, SIGNAL(mapped(const Wt::WString &)),this, SLOT(SaveIconCallback(const Wt::WString&)));
    Wt::WString txt = aCommand + fStringSeparator + aLabel;
    signalMapper->setMapping(action, Wt::WString(txt));
    
    // special cases : cursor style
  } else if ((std::string(aIconFile) == "move") ||
             (std::string(aIconFile) == "rotate") ||
             (std::string(aIconFile) == "pick") ||
             (std::string(aIconFile) == "zoom_out") ||
             (std::string(aIconFile) == "zoom_in")) {
    action->setCheckable(TRUE);
    action->setChecked(TRUE);
    action->setData(aIconFile);
    
    connect(signalMapper, SIGNAL(mapped(const Wt::WString &)),this, SLOT(ChangeCursorStyle(const Wt::WString&)));
    signalMapper->setMapping(action, Wt::WString(aIconFile));
    
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
    connect(signalMapper, SIGNAL(mapped(const Wt::WString &)),this, SLOT(ChangeSurfaceStyle(const Wt::WString&)));
    signalMapper->setMapping(action, Wt::WString(aIconFile));
    
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
    connect(signalMapper, SIGNAL(mapped(const Wt::WString &)),this, SLOT(ChangePerspectiveOrthoCallback(const Wt::WString&)));
    signalMapper->setMapping(action, Wt::WString(aIconFile));
    
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
    
    connect(signalMapper, SIGNAL(mapped(const char *)),this, SLOT(ButtonCallback(const char *)));
    signalMapper->setMapping(action, aCommand);
  }
*/
}


/** Create a widget with the command parameters inside
 @param command: command line
 @parent : parent widget
 @isDialog : true if we want apply/cancel button and close at end, false if we want only apply
 */
bool G4UIWt::CreateCommandWidget(G4UIcommand* aCommand, Wt::WContainerWidget* aParent, bool isDialog) {
  
  if (aCommand == NULL) {
    return false;
  }
  
  
  // parameters
  G4int n_parameterEntry = aCommand->GetParameterEntries();
  if( n_parameterEntry > 0 ) {
    G4UIparameter *param;
    
    // Re-implementation of G4UIparameter.cc
    Wt::WContainerWidget* paramWidget = new Wt::WContainerWidget();
    Wt::WGridLayout* gridLayout = new Wt::WGridLayout(paramWidget);
    
    // Special case for colour, try to display a color chooser if we found red/green/blue parameter
    unsigned int nbColorParameter = 0;
    bool isStillColorParameter = false;
    bool isColorDialogAdded = false;
    Wt::WLabel* redLabel = NULL;
    Wt::WLabel* greenLabel = NULL;
    Wt::WString redDefaultStr = "";
    Wt::WString greenDefaultStr = "";
    Wt::WString blueDefaultStr = "";
    Wt::WWidget* redInput = NULL;
    Wt::WWidget* greenInput = NULL;
    
    for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ ) {
      Wt::WString txt;
      param = aCommand->GetParameter(i_thParameter);
      Wt::WLabel* label = new Wt::WLabel(Wt::WString((char*)(param->GetParameterName()).data()));
      
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
      Wt::WWidget* input = NULL;
      char paramType = param->GetParameterType();
      if ((paramType == 'd') || (paramType == 'i')) {
        input = new Wt::WLineEdit();
        // set default value
        dynamic_cast<Wt::WLineEdit*>(input)->setText(Wt::WString((char*)(param->GetDefaultValue()).data()));
        
        if (((label->text() == "red") || (label->text() == "red_or_string")) && isStillColorParameter) {
          redDefaultStr = Wt::WString((char*)(param->GetDefaultValue()).data());
        } else if ((label->text() == "green") && isStillColorParameter) {
          greenDefaultStr = Wt::WString((char*)(param->GetDefaultValue()).data());
        } else if ((label->text() == "green") && isStillColorParameter) {
          blueDefaultStr = Wt::WString((char*)(param->GetDefaultValue()).data());
        }
        
      } else if (paramType == 'b') {
        input = new Wt::WContainerWidget();
        Wt::WHBoxLayout* layout = new Wt::WHBoxLayout(input);
        
        Wt::WButtonGroup* buttons = new Wt::WButtonGroup();
        Wt::WRadioButton* radioOff = new Wt::WRadioButton("0");
        Wt::WRadioButton* radioOn = new Wt::WRadioButton("1");
        buttons->addButton(radioOn);
        buttons->addButton(radioOff);
        layout->addWidget(radioOn);
        layout->addWidget(radioOff);
        
        // set default value
        Wt::WString defaultValue = Wt::WString((char*)(param->GetDefaultValue()).data());
        if (defaultValue == "0") {
          radioOff->setChecked(true);
        } else if (defaultValue == "1") {
          radioOn->setChecked(true);
        }
      } else if ((paramType == 's') && (!param->GetParameterCandidates().isNull())) {
        input = new Wt::WComboBox();
        Wt::WString candidates = Wt::WString((char*)(param->GetParameterCandidates()).data());
        printf("*** G4UIWt::CreateCommandWidget, missing parameter management for 's'\n");
        /*        QStringList list = candidates.split (" ");
        
        // add all candidates to widget
        Wt::WString defaultValue = Wt::WString((char*)(param->GetDefaultValue()).data());
        for (int a=0; a<list.size(); a++) {
          dynamic_cast<Wt::WComboBox*>(input)->addItem(list.at(a));
          if (list.at(a) == defaultValue) {
            dynamic_cast<Wt::WComboBox*>(input)->setCurrentIndex(a);
          }
        }
*/ 
      } else if (paramType == 's') {  // string
        input = new Wt::WLineEdit();
        // set default value
        dynamic_cast<Wt::WLineEdit*>(input)->setText(Wt::WString((char*)(param->GetDefaultValue()).data()));
        
      } else if (paramType == 'c') {  // on/off
        input = new Wt::WContainerWidget();
        Wt::WHBoxLayout* layout = new Wt::WHBoxLayout(input);
        
        Wt::WButtonGroup* buttons = new Wt::WButtonGroup();
        Wt::WRadioButton* radioOff = new Wt::WRadioButton("off");
        Wt::WRadioButton* radioOn = new Wt::WRadioButton("on");
        buttons->addButton(radioOn);
        buttons->addButton(radioOff);
        layout->addWidget(radioOn);
        layout->addWidget(radioOff);
        
        // set default value
        Wt::WString defaultValue = Wt::WString((char*)(param->GetDefaultValue()).data());
        if (defaultValue == "off") {
          radioOff->setChecked(true);
        } else if (defaultValue == "on") {
          radioOn->setChecked(true);
        }
        
      } else {
        input = new Wt::WLineEdit();
        dynamic_cast<Wt::WLineEdit*>(input)->setText(Wt::WString((char*)(param->GetDefaultValue()).data()));
      }
      
      txt += "\nParameter : " + Wt::WString((char*)(param->GetParameterName()).data()) + "\n";
      if( ! param->GetParameterGuidance().isNull() )
        txt += Wt::WString((char*)(param->GetParameterGuidance()).data())+ "\n" ;
      
// FIXME ?
      txt += Wt::WString(" Parameter type  : ") + std::string(&paramType).c_str() + "\n";
      if(param->IsOmittable()){
        txt += " Omittable       : True\n";
      } else {
        txt += " Omittable       : False\n";
      }
      if( param->GetCurrentAsDefault() ) {
        txt += " Default value   : taken from the current value\n";
      } else if( ! param->GetDefaultValue().isNull() ) {
        txt += " Default value   : " + Wt::WString((char*)(param->GetDefaultValue()).data())+ "\n";
      }
      if( ! param->GetParameterRange().isNull() ) {
        txt += " Parameter range : " + Wt::WString((char*)(param->GetParameterRange()).data())+ "\n";
      }
      if( ! param->GetParameterCandidates().isNull() ) {
        txt += " Candidates      : " + Wt::WString((char*)(param->GetParameterCandidates()).data())+ "\n";
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
          Wt::WColor wc;
          if ((redDefaultStr != "") && (redDefaultStr != "") && (redDefaultStr != "")) {
            // 255 max
            wc.setRgb(atof(redDefaultStr.toUTF8().c_str())*256,
                       atof(greenDefaultStr.toUTF8().c_str())*256,
                       atof(blueDefaultStr.toUTF8().c_str())*256);
          }
          printf("*** G4UIWt::CreateCommandWidget, missing icon on command widget\n");
/*
 QPixmap pixmap = QPixmap(QSize(16, 16));
          pixmap.fill (wc);
          Wt::WPainter painter(&pixmap);
          painter.setPen(Wt::black);
          painter.drawRect(0,0,15,15); // Draw contour
          
          input = new Wt::WPushButton("Change color");
          dynamic_cast<Wt::WPushButton*>(input)->setIcon(pixmap);
          dynamic_cast<Wt::WPushButton*>(input)->setAccessibleName(redDefaultStr+" "+greenDefaultStr+" "+blueDefaultStr);
          label = new Wt::WLabel("Choose color");
          
          // less 1 because we have to add one to the row number
          nbColorParameter--;
          gridLayout->addWidget(label,i_thParameter-nbColorParameter,0);
          input->setToolTip("Select the current color");
          gridLayout->addWidget(input,i_thParameter-nbColorParameter,1);
          
          // Connect pushButton to ColorDialog in callback
          QSignalMapper* signalMapper = new QSignalMapper(this);
          signalMapper->setMapping(input,input);
          connect(input, SIGNAL(clicked()), signalMapper, SLOT(map()));
          connect(signalMapper, SIGNAL(mapped(Wt::WWidget*)),this, SLOT(ChangeColorCallback(Wt::WWidget*)));
          
*/
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
    Wt::WLabel* name = new Wt::WLabel(Wt::WString((char*)(aCommand->GetCommandPath().data())));
    name->hide();
    gridLayout->addWidget(name,n_parameterEntry-nbColorParameter,0);
    
    Wt::WPushButton* applyButton = new Wt::WPushButton("Apply");
    if (!isDialog) {
      
      gridLayout->addWidget(applyButton,n_parameterEntry-nbColorParameter,1);
      
      printf("*** G4UIWt::CreateCommandWidget, missing connection on ApplyButton\n");
/*
 QSignalMapper* signalMapper = new QSignalMapper(this);
      signalMapper->setMapping(applyButton, paramWidget);
      connect(applyButton, SIGNAL(clicked()), signalMapper, SLOT(map()));
      connect(signalMapper, SIGNAL(mapped(Wt::WWidget*)),this, SLOT(VisParameterCallback(Wt::WWidget*)));
*/
    } else {
      // Apply/Cancel buttons
      
      printf("*** G4UIWt::CreateCommandWidget, missing connection on Apply/Cancel Button\n");
/*      applyButton->setAutoDefault( TRUE );
      applyButton->setDefault( TRUE );
      gridLayout->addWidget(applyButton,n_parameterEntry-nbColorParameter,0);
      
      Wt::WPushButton* cancelButton = new Wt::WPushButton( tr( "&Cancel" ));
      cancelButton->setAutoDefault( TRUE );
      gridLayout->addWidget(cancelButton,n_parameterEntry-nbColorParameter,1);
      
      QSignalMapper* signalMapper = new QSignalMapper(this);
      signalMapper->setMapping(applyButton, paramWidget);
      connect(applyButton, SIGNAL(clicked()), signalMapper, SLOT(map()));
      connect(signalMapper, SIGNAL(mapped(Wt::WWidget*)),this, SLOT(VisParameterCallback(Wt::WWidget*)));
      
      Wt::WWidget * parentCheck = aParent;
      Wt::WDialog* parentDialog = NULL;
      bool found = false;
      while ((parentCheck->parentWidget()) != NULL) {
        parentCheck = (Wt::WWidget*) parentCheck->parentWidget();
        parentDialog = dynamic_cast<Wt::WDialog*>(parentCheck);
        if (parentDialog) {
          connect( applyButton, SIGNAL( clicked() ), parentDialog, SLOT( accept() ) );
          connect( cancelButton, SIGNAL( clicked() ), parentDialog, SLOT( reject() ) );
          found = true;
        }
      }
      if (!found) {
        return false;
      }
*/
    }
    
    if (!aParent->layout()) {
      aParent->setLayout(new Wt::WVBoxLayout());
    }
    aParent->layout()->addWidget(paramWidget);
  }
  
  return true;
}


/**   Event filter method. Every event from WtApplication goes here.<br/>
 We apply a filter only for the Up and Down Arrow press when the Wt::WLineEdit<br/>
 is active. If this filter match, Up arrow we give the previous command<br/>
 and Down arrow will give the next if exist.<br/>
 @param obj Emitter of the event
 @param event Kind of event
 */
bool G4UIWt::eventFilter( // Should stay with a minuscule eventFilter because of Wt
                         Wt::WObject* /* aObj */
                         ,Wt::WEvent* /* aEvent */
                         )
{
  bool res= false;
  printf("*** G4UIWt::eventFilter, missing eventFilter on everything\n");
/*
 bool moveCommandCursor = false;
  if (aObj == NULL) return false;
  if (aEvent == NULL) return false;
  
  if (aObj == fHistoryTBTableList) {
    if (aEvent->type() == Wt::WEvent::KeyPress) {
      fCommandArea->setFocus();
    }
  }
  if (aObj == fCommandArea) {
    if (aEvent->type() == Wt::WEvent::KeyPress) {
      QKeyEvent *e = static_cast<QKeyEvent*>((QEvent*)aEvent);
      if ((e->key() == (Wt::Key_Down)) ||
          (e->key() == (Wt::Key_PageDown)) ||
          (e->key() == (Wt::Key_Up)) ||
          (e->key() == (Wt::Key_PageUp))) {
        int selection = fHistoryTBTableList->currentRow();
        if (fHistoryTBTableList->count()) {
          if (selection == -1) {
            selection = fHistoryTBTableList->count()-1;
          } else {
            if (e->key() == (Wt::Key_Down)) {
              if (selection <(fHistoryTBTableList->count()-1))
                selection++;
            } else if (e->key() == (Wt::Key_PageDown)) {
              selection = fHistoryTBTableList->count()-1;
            } else if (e->key() == (Wt::Key_Up)) {
              if (selection >0)
                selection --;
            } else if (e->key() == (Wt::Key_PageUp)) {
              selection = 0;
            }
          }
          fHistoryTBTableList->clearSelection();
          fHistoryTBTableList->item(selection)->setSelected(true);
          fHistoryTBTableList->setCurrentItem(fHistoryTBTableList->item(selection));
        }
        moveCommandCursor = true;
      } else if (e->key() == (Wt::Key_Tab)) {
        G4String ss = Complete(fCommandArea->text().toUTF8().c_str());
        fCommandArea->setText((char*)(ss.data()));
        
        // do not pass by parent, it will disable widget tab focus !
        return true;
        // L.Garnier : MetaModifier is CTRL for MAC, but I don't want to put a MAC
        // specific #ifdef
      } else if (((e->modifiers () == Wt::ControlModifier) || (e->modifiers () == Wt::MetaModifier)) && (e->key() == Wt::Key_A)) {
        fCommandArea->home(false);
        return true;
      } else if (((e->modifiers () == Wt::ControlModifier) || (e->modifiers () == Wt::MetaModifier)) && (e->key() == Wt::Key_E)) {
        fCommandArea->end(false);
        return true;
      }
    }
  }
  // change cursor position if needed
  if (moveCommandCursor == true) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::eventFilter setCursor Position\n");
#endif
    fCommandArea->setCursorPosition ( fCommandArea->text().length() );
    fCommandArea->setCursorPosition (4);
  } else {
    // pass the event on to the parent class
    res = QObject::eventFilter(aObj, aEvent);
  }
*/
  return res;
}


/**   This callback is activated when user selected a item in the help tree
 */
void G4UIWt::HelpTreeClicCallback (
)
{
  Wt::WTreeNode* item =  NULL;
  if (!fHelpTreeWidget)
    return ;
  
  if (!fHelpArea)
    return;
  
  const Wt::WTree::WTreeNodeSet& list = fHelpTreeWidget->selectedNodes();
  if (list.empty())
    return;
  item = *list.begin();
  if (!item)
    return;
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();
  
  std::string itemText = GetLongCommandPath(item).toUTF8();
  
  G4UIcommand* command = treeTop->FindPath(itemText.c_str());
  
  if (command) {
    fHelpArea->setText(GetCommandList(command));
  } else {  // this is a command
    G4UIcommandTree* path = treeTop->FindCommandTree(itemText.c_str());
    if ( path) {
      // this is not a command, this is a sub directory
      // We display the Title
      fHelpArea->setText(path->GetTitle().data());
    }
  }
}

/**   This callback is activated when user double clic on a item in the help tree
 */
void G4UIWt::HelpTreeDoubleClicCallback (
)
{
  HelpTreeClicCallback();
  
  Wt::WTreeNode* item =  NULL;
  if (!fHelpTreeWidget)
    return ;
  
  if (!fHelpArea)
    return;
  
  const Wt::WTree::WTreeNodeSet& list = fHelpTreeWidget->selectedNodes();
  if (list.empty())
    return;
  item = *list.begin();
  if (!item)
    return;
  
  fCommandArea->setText("");
  fCommandArea->setText(GetLongCommandPath(item));
}


/**
 Receive a cout from Geant4. We have to display it in the cout zone
 @param aString : label to add in the display area
 @return 0
 */
G4int G4UIWt::ReceiveG4cout (
                             const G4String& aString
                             )
{
  if (!aString) return 0;
  
  Wt::WStringListModel newStr;
  
  // Add to stringList
  std::string whiteSpaces( " \f\n\r\t\v" );
  std::string path = (char*)aString.data();
  
  std::string::size_type posR = path.find_last_not_of( whiteSpaces );
  path.erase( posR + 1 );
  
  std::string::size_type posL = path.find_first_not_of( whiteSpaces );
  path.erase( 0, posL );

  printf("*** G4UIWt::ReceiveG4cout, missing filtering\n");
/*  newStr = Wt::WString(path);
  fG4cout += newStr;
  
  QStringList result = newStr.filter(fCoutFilter->text());
 
  if (result.join("").isEmpty()) {
    return 0;
  }
 */
  if (!fCoutTBTextArea) {
    printf("*** G4UIWt::ReceiveG4cout, create a new fCoutTBTextArea \n");
    fCoutTBTextArea = new Wt::WTextArea();
    fCoutTBTextArea->setText("");
  }
  fCoutTBTextArea->setText(fCoutTBTextArea->text()+"\n"+path);
  fCoutTBTextArea->refresh();
  
  return 0;
}


/**
 Receive a cerr from Geant4. We have to display it in the cout zone
 @param aString : label to add in the display area
 @return 0
 */
G4int G4UIWt::ReceiveG4cerr (
                             const G4String& aString
                             )
{
  if (!aString) return 0;
  
  Wt::WStringListModel newStr;
  
  // Add to stringList
  std::string whiteSpaces( " \f\n\r\t\v" );
  std::string path = (char*)aString.data();
  
  std::string::size_type posR = path.find_last_not_of( whiteSpaces );
  path.erase( posR + 1 );
  
  std::string::size_type posL = path.find_first_not_of( whiteSpaces );
  path.erase( 0, posL );
  
  printf("*** G4UIWt::ReceiveG4cerr, missing filtering\n");
  /*   newStr = Wt::WStringList(path);
  fG4cout += newStr;
  
   QStringList result = newStr.filter(fCoutFilter->text());
  
   */
  // Suppress space, \n,\t,\r...
  if (path != "") {
    if ((G4StateManager::GetStateManager()->GetCurrentState() == G4State_Abort) ||
        (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Quit )) {
      // In case of Abort or Quit, the useful error message should be in the last error message !
      Wt::WMessageBox::show("Error", Wt::WString(fLastErrMessage.data())+"\n"+aString.data(), Wt::Ok );
    }
  
  }
  printf("****ERR****: %s\n",path.c_str());
  if (!fCoutTBTextArea) {
    fCoutTBTextArea = new Wt::WTextArea();
  }
  fCoutTBTextArea->setText(fCoutTBTextArea->text()+"\n<font color='red'>"+path+"</font>");
  fCoutTBTextArea->refresh();

  
  if (path != "") {
    fLastErrMessage = aString;
  }
  return 0;
}


/** Callback when the text in the line edit is changed.
 When a newline is inserted, trigger the Activate Command
 on this text end set unchanged the end of the line after the newline.
 */
void G4UIWt::CommandEditedCallback(const Wt::WString &)
{
  printf("*** G4UIWt::CommandEditedCallback, missing callback on command line edit\n");
/*  QStringList list = fCommandArea->text().split(QRegExp("[\r\n]"),QString::SkipEmptyParts);
  
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
*/
}



void G4UIWt::CoutFilterCallback(
                                const Wt::WString & /* text */) {
  printf("*** G4UIWt::CoutFilterCallbackt, missing filtering\n");
  
/*  QStringList result = fG4cout.filter(text);
  fCoutTBTextArea->setPlainText(result.join("\n"));
  
  fCoutTBTextArea->repaint();
  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
  
*/
}





/**
 Add a new menu to the menu bar
 @param aName name of menu
 @param aLabel label to display
 */
void G4UIWt::AddMenu (
                      const char* aName
                      ,const char* aLabel
                      )
{
  if (aName == NULL) return;
  if (aLabel == NULL) return;
  
  printf("*** G4UIWt::AddMenu, missing \n");
/*  Wt::WMenu *fileMenu = new Wt::WMenu(aLabel);
  fMainWindow->menuBar()->addMenu(fileMenu);
  
  AddInteractor (aName,(G4Interactor)fileMenu);
*/
}


/**
 Add the following command to the corresponding groupbox
 If depthLevel is 1 : create ToolBox
 If depthLevel is 2 or more : create GroupBox
 */
bool G4UIWt::CreateVisCommandGroupAndToolBox(
                                             G4UIcommand* /* aCommand */
                                             ,Wt::WWidget* /* aParent */
                                             ,int /* aDepthLevel */
                                             ,bool /* isDialog */
                                             )
{
  printf("*** G4UIWt::CreateVisCommandGroupAndToolBox, missing \n");
  /*
   std::string str ((char*)(aCommand->GetCommandPath().data()));
  std::string str2 ("/");
  std::size_t pos = -1;
  
  for (int a=1; a<=-aDepthLevel; a++) {
    pos = str.find(str2, pos+1);
  }
  
  Wt::WString commandText2;
  std::string commandText = "";
  if (pos!=std::string::npos) {
    commandText = str.substr(pos);
  }
  // FIXME :  commandText2 = Wt::WString((char*)(aCommand->GetCommandPath().data())).section("/",-aDepthLevel);
  // FIXME :   printf(" CommandText : %s\n",commandText.toUTF8().c_str());
  // FIXME :   printf(" CommandText2 : %s\n",commandText2.toUTF8().c_str());
  
  if (commandText == "") {
    return false;
  }
  
  // Look if groupBox is create
  //  Wt::WGroupBox* gBoxCommandWidget;
  Wt::WWidget* newParentWidget = NULL;
  bool found = false;
  
  pos = commandText.find("/", 0);
  Wt::WString commandSection = "";
  if (pos!=std::string::npos) {
    commandSection = Wt::WString(str.substr(0,pos).c_str());
  }
  
  // FIXME :   Wt::WString commandSection = commandText.left(commandText.indexOf("/"));
  
  if (aDepthLevel == 1) {
    Wt::WToolBox* currentParent = dynamic_cast<Wt::WToolBox*>(aParent);
    if (currentParent != 0){
      
      // already exists ?
      for (int a=0; a<currentParent->count(); a++) {
        if (currentParent->itemText(a) == commandSection) {
          found = true;
          newParentWidget = (Wt::WWidget*) currentParent->widget(a);
        }
      }
    }
    // Not found ? create it
    if (!found) {
      newParentWidget = new Wt::WGroupBox();
      //        newParentWidget->setSizePolicy (QSizePolicy(QSizePolicy::Maximum,QSizePolicy::Maximum));
      newParentWidget->setLayout(new Wt::WVBoxLayout());
      if (currentParent != 0){
        currentParent->addItem(newParentWidget,commandSection);
      } else {
        if (!aParent->layout()) {
          aParent->setLayout(new Wt::WVBoxLayout());
        }
        aParent->layout()->addWidget(newParentWidget);
      }
      
      if (commandText.find("/", 0) != std::string::npos) {
        
        // Guidance
        Wt::WString guidance;
        G4int n_guidanceEntry = aCommand->GetGuidanceEntries();
        for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
          guidance += Wt::WString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
        }
        newParentWidget->setToolTip(guidance);
      }
      
      Wt::WScrollArea* sc = dynamic_cast<Wt::WScrollArea*>(newParentWidget->parent()->parent());
      if (sc != 0) {
        sc->ensureWidgetVisible(newParentWidget);
        //          sc->setSizePolicy (QSizePolicy(QSizePolicy::Maximum,QSizePolicy::Maximum));
        
      }
    }
  } else {
    
    // try to know if this level is already there
    Wt::WGroupBox* currentParent = dynamic_cast<Wt::WGroupBox*>(aParent);
    if (currentParent != 0){
      
      // if depth==2, then we add a [more parameters inside] to the toolBoxItem parent
      // Wt::WGroupBox > Wt::WWidget > Wt::WScrollArea > Wt::WToolBox
      if (aDepthLevel == 2){
        Wt::WToolBox* parentToolBox = dynamic_cast<Wt::WToolBox*>(((Wt::WWidget*) (currentParent))->parent()->parent()->parent());
        if (parentToolBox != 0) {
          //          parentToolBox->setItemText(parentToolBox->indexOf(currentParent),"[more parameters inside]");
        }
      }
      for (int a=0; a<aParent->layout()->count(); a++) {
        Wt::WGroupBox* gb = dynamic_cast<Wt::WGroupBox*>(aParent->layout()->itemAt(a)->widget());
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
      newParentWidget = new Wt::WGroupBox(commandSection);
      newParentWidget->setLayout(new Wt::WVBoxLayout());
      if (!aParent->layout()) {
        aParent->setLayout(new Wt::WVBoxLayout());
      }
      aParent->layout()->addWidget(newParentWidget);
      
      // set toolTip
      // Guidance
      Wt::WString guidance;
      G4int n_guidanceEntry = aCommand->GetGuidanceEntries();
      for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
        guidance += Wt::WString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
      }
      newParentWidget->setToolTip(guidance);
    }
  }
  
  // fill command groupbox
  if (commandText.find("/", 0) != std::string::npos) {
    if (CreateCommandWidget(aCommand, newParentWidget,isDialog)) {
      return true;
    }
  } else {
    CreateVisCommandGroupAndToolBox(aCommand,newParentWidget, aDepthLevel-1,isDialog);
  }
*/  
  return true;
}


/** Callback when one of the scene/vis parameters has changed
 */
void G4UIWt::VisParameterCallback(Wt::WContainerWidget* widget){
  if (widget == NULL) {
    return;
  }
  printf("*** G4UIWt::VisParameterCallback, missing \n");
/*
  // Look in all the Grid layout, but only column 1 (0 is the parameter name)
  Wt::WGridLayout* grid = dynamic_cast<Wt::WGridLayout*>(widget->layout());
  if (grid == 0) {
    return;
  }
  Wt::WString command;
  Wt::WWidget* name = (Wt::WWidget*) (grid->itemAtPosition(grid->rowCount()-1,0)->widget());
  if (widget == NULL) {
    return;
  }
  if (dynamic_cast<Wt::WLabel*>(name) == 0) {
    return;
  }
  command += (dynamic_cast<Wt::WLabel*>(name))->text()+" ";
  
  for (int a=0;a<grid->rowCount()-1; a++) {
    Wt::WWidget* widgetTmp = (Wt::WWidget*) (grid->itemAtPosition(a,1)->widget());
    
    // 4 kind of widgets : Wt::WLineEdit / Wt::WComboBox / radioButtonsGroup / Wt::WPushButton (color chooser)
    if (widgetTmp != NULL) {
      
      if (dynamic_cast<Wt::WLineEdit*>(widgetTmp) != 0) {
        command += (dynamic_cast<Wt::WLineEdit*>(widgetTmp))->text()+" ";
        
      } else if (dynamic_cast<Wt::WComboBox*>(widgetTmp) != 0){
        command += (dynamic_cast<Wt::WComboBox*>(widgetTmp))->itemText((dynamic_cast<Wt::WComboBox*>(widgetTmp))->currentIndex())+" ";
        
        // Color chooser
      } else if (dynamic_cast<Wt::WPushButton*>(widgetTmp) != 0){
        command += widgetTmp->accessibleName()+" ";
        
        // Check for Button group
      } else if (dynamic_cast<Wt::WWidget*>(widgetTmp) != 0){
        if (widgetTmp->layout()->count() > 0){
          if (dynamic_cast<Wt::WRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()) != 0) {
            QAbstractButton * checked = (dynamic_cast<Wt::WRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()))->group()->checkedButton();
            if (checked != 0) {
              command += (dynamic_cast<Wt::WRadioButton*>(widgetTmp->layout()->itemAt(0)->widget()))->group()->checkedButton()->text()+" ";
            }
          }
        }
        
      }
    }
  }
  if (command != "") {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if(UI != NULL)  {
      UI->ApplyCommand(command.toUTF8().c_str());
    }
  }
 */
}


/**   Callback called when user select an old command in the command history<br>
 Give it to the command area.
 */
void G4UIWt::CommandHistoryCallback(
)
{
  if (!fHistoryTBTableList)
    return ;
  fCommandArea->setText(fHistoryTBTableList->currentText ());
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CommandHistoryCallback change text\n");
#endif
}


void G4UIWt::OpenHelpTreeOnCommand(
                                   const Wt::WString & /* searchText */
                                   )
{
  printf("*** G4UIWt::OpenHelpTreeOnCommand, missing \n");
/*  // the help tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();
  
  G4int treeSize = treeTop->GetTreeEntry();
  
  // clear old help tree
  fHelpTreeWidget->clear();
  
  // look for new items
  
  int tmp = 0;
  
  QMap<int,Wt::WString> commandResultMap;
  QMap<int,Wt::WString> commandChildResultMap;
  
  for (int a=0;a<treeSize;a++) {
    G4UIcommand* command = treeTop->FindPath(treeTop->GetTree(a+1)->GetPathName().data());
    tmp = GetCommandList (command).count(searchText,Wt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,Wt::WString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(treeTop->GetTree(a+1),searchText);
    // insert new childs
    if (!commandChildResultMap.empty()) {
      QMap<int,Wt::WString>::const_iterator i = commandChildResultMap.constBegin();
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
  Wt::WStringList labels;
  labels << Wt::WString("Command") << Wt::WString("Match");
  fHelpTreeWidget->setHeaderLabels(labels);
  
  if (commandResultMap.empty()) {
    fHelpArea->setText("No match found");
    return;
  }
  
  QMap<int,Wt::WString>::const_iterator i = commandResultMap.constEnd();
  i--;
  // 10 maximum progress values
  float multValue = 10.0/(float)(i.key());
  Wt::WString progressChar = "|";
  Wt::WString progressStr = "|";
  
  Wt::WTreeNode * newItem;
  bool end = false;
  while (!end) {
    if (i == commandResultMap.constBegin()) {
      end = true;
    }
    for(int a=0;a<int(i.key()*multValue);a++) {
      progressStr += progressChar;
    }
    newItem = new Wt::WTreeNode(fHelpTreeWidget);
    Wt::WString commandStr = i.value().trimmed();
    
    if (commandStr.indexOf("/") == 0) {
      commandStr = commandStr.right(commandStr.size()-1);
    }
    
    newItem->setText(0,commandStr);
    newItem->setText(1,progressStr);
    
    newItem->setForeground ( 1, QBrush(Wt::blue) );
    progressStr = "|";
    i--;
  }
  fHelpTreeWidget->resizeColumnToContents (0);
  fHelpTreeWidget->sortItems(1,Wt::DescendingOrder);
  //  fHelpTreeWidget->setColumnWidth(1,10);//resizeColumnToContents (1);
*/
}

/*
WMap<int,Wt::WString> G4UIWt::LookForHelpStringInChildTree(
                                                                   G4UIcommandTree *aCommandTree
                                                                   ,const Wt::WString & text
                                                                   )
{
  QMap<int,Wt::WString> commandResultMap;
  if (aCommandTree == NULL) return commandResultMap;
  
  
  // Get the Sub directories
  int tmp = 0;
  QMap<int,Wt::WString> commandChildResultMap;
  
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {
    const G4UIcommand* command = aCommandTree->GetGuidance();
    tmp = GetCommandList (command).count(text,Wt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,Wt::WString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()));
    }
    // look for childs
    commandChildResultMap = LookForHelpStringInChildTree(aCommandTree->GetTree(a+1),text);
    
    if (!commandChildResultMap.empty()) {
      // insert new childs
      QMap<int,Wt::WString>::const_iterator i = commandChildResultMap.constBegin();
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
    tmp = GetCommandList (command).count(text,Wt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,Wt::WString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()));
    }
    
  }
  return commandResultMap;
}
*/


void G4UIWt::ChangeColorCallback(Wt::WContainerWidget* widget) {
  if (widget == NULL) {
    return;
  }
  printf("*** G4UIWt::ChangeColorCallback, missing \n");
/*
  Wt::WPushButton* button = dynamic_cast<Wt::WPushButton*>(widget);
  if (button == 0) {
    return;
  }
  Wt::WString value = button->accessibleName();
  
  Wt::WColor old;
  old.setRgbF(value.section(" ",0,1).toDouble(),
              value.section(" ",1,2).toDouble(),
              value.section(" ",2,3).toDouble());
  Wt::WColor color = Wt::WColorDialog::getColor(old,
                                                    fUITabWidget,
                                                    "Change color",
                                                    Wt::WColorDialog::ShowAlphaChannel);
  
  
  if (color.isValid()) {
    // rebuild the widget icon
    QPixmap pixmap = QPixmap(QSize(16, 16));
    pixmap.fill (color);
    Wt::WPainter painter(&pixmap);
    painter.setPen(Wt::black);
    painter.drawRect(0,0,15,15); // Draw contour
    
    button->setAccessibleName(Wt::WString::number(color.redF())+" "+
                              Wt::WString::number(color.greenF())+" "+
                              Wt::WString::number(color.blueF())+" "
                              );
    button->setIcon(pixmap);
    
    
  }
*/
}

void G4UIWt::ChangeCursorStyle(const Wt::WString& /* action */) {
  
  // Theses actions should be in the app toolbar
  
  fMoveSelected = true;
  fPickSelected = true;
  fRotateSelected = true;
  fZoomInSelected = true;
  fZoomOutSelected = true;
  
  printf("*** G4UIWt::ChangeCursorStyle, missing \n");
/*
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
*/
}


/* A little bit like "void G4OpenGLWtViewer::toggleDrawingAction(int aAction)"
 But for all viewers, not only Wt
 
 FIXME : Should be a feedback when changing viewer !
 
 */
void G4UIWt::ChangeSurfaceStyle(const Wt::WString& /* action */) {
  
  // Theses actions should be in the app toolbar
  
  printf("*** G4UIWt::ChangeSurfaceStyle, missing \n");
/*  if (fToolbarApp == NULL) return;
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
*/
}

void G4UIWt::OpenIconCallback(const Wt::WString& /* aParam */) {
  
  printf("*** G4UIWt::OpenIconCallback, missing \n");
/*
 Wt::WString aCommand = aParam.left(aParam.indexOf(fStringSeparator));
  Wt::WString aLabel = aParam.mid(aParam.indexOf(fStringSeparator)+fStringSeparator.length());
  
  Wt::WString nomFich = Wt::WFileDialog::getOpenFileName(fMainWindow, aLabel, fLastOpenPath, "Macro files (*.mac)");
  if (nomFich != "") {
    G4UImanager::GetUIpointer()->ApplyCommand((Wt::WString(aCommand)+ Wt::WString(" ")+ nomFich).toUTF8().c_str());
    QDir dir;
    fLastOpenPath = dir.absoluteFilePath(nomFich);
  }
*/
}


void G4UIWt::SaveIconCallback(const Wt::WString& /* aParam */) {
  
  printf("*** G4UIWt::SaveIconCallback, missing \n");

  /*
   Wt::WString aCommand = aParam.left(aParam.indexOf(fStringSeparator));
  Wt::WString aLabel = aParam.mid(aParam.indexOf(fStringSeparator)+fStringSeparator.length());
  
  Wt::WString nomFich = Wt::WFileDialog::getSaveFileName(fMainWindow, aLabel, fLastOpenPath, "Macro files (*.mac)");
  if (nomFich != "") {
    G4UImanager::GetUIpointer()->ApplyCommand((Wt::WString(aCommand)+ Wt::WString(" ")+nomFich).toUTF8().c_str());
    QDir dir;
    fLastOpenPath = dir.absoluteFilePath(nomFich);
  }
*/
}


void G4UIWt::ChangePerspectiveOrthoCallback(const Wt::WString& /* action */) {
  
  // Theses actions should be in the app toolbar
  
  printf("*** G4UIWt::ChangePerspectiveOrthoCallback, missing \n");
/*
 if (fToolbarApp == NULL) return;
  QList<QAction *> list = fToolbarApp->actions ();
  Wt::WString checked = "";
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
*/
}

/*
void G4UIWt::SetIconMoveSelected() {
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


void G4UIWt::SetIconRotateSelected() {
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


void G4UIWt::SetIconPickSelected() {
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


void G4UIWt::SetIconZoomInSelected() {
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


void G4UIWt::SetIconZoomOutSelected() {
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


void G4UIWt::SetIconSolidSelected() {
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


void G4UIWt::SetIconWireframeSelected() {
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


void G4UIWt::SetIconHLRSelected() {
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


void G4UIWt::SetIconHLHSRSelected() {
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


void G4UIWt::SetIconPerspectiveSelected() {
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



void G4UIWt::SetIconOrthoSelected() {
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
*/


/**
 Create the help tree widget
 @param parent : parent of tree widget
 @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIWt::FillHelpTree()
{
  if (! fHelpTreeWidget ) {
    InitHelpTreeAndVisParametersWidget();
  }
  
  Wt::WString searchText = fHelpLine->text();
  
  if (searchText =="") {
    // clear old help tree
    //    fHelpTreeWidget->clear();
  } else {
    return;
  }
  
  if (fHelpArea) {
    fHelpArea->setText("");
  }
  
  if (fHelpLine) {
    fHelpLine->setText("");
  }
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();
  
  G4int treeSize = treeTop->GetTreeEntry();
  Wt::WTreeNode * newItem = NULL;
  Wt::WString commandText = "";
  
  for (int a=0;a<treeSize;a++) {
    // Creating new item
    newItem = NULL;
    
    // trim path
    std::string whiteSpaces( " \f\n\r\t\v" );
    std::string path = (treeTop->GetTree(a+1)->GetPathName()).data();
    
    std::string::size_type posR = path.find_last_not_of( whiteSpaces );
    path.erase( posR + 1 );
    
    std::string::size_type posL = path.find_first_not_of( whiteSpaces );
    path.erase( 0, posL );
    
    commandText = Wt::WString(path.c_str());
    
    // if already exist, don't create it !
    if (fHelpTreeWidget->treeRoot()) {
      for (int b=0;b<fHelpTreeWidget->treeRoot()->displayedChildCount();b++) {
      if (!newItem)
          newItem = FindTreeItem(fHelpTreeWidget->treeRoot()->childNodes()[b],path.c_str());
      }
    }
    if (newItem == NULL) {
      newItem = new Wt::WTreeNode(GetShortCommandPath(path), 0);
      fHelpTreeWidget->setTreeRoot(newItem);
    }
    
    // look for childs
    CreateHelpTree(newItem,treeTop->GetTree(a+1));
  }
  
}



/**
 Called by intercoms/src/G4UImanager.cc<br>
 Called by visualization/management/src/G4VisCommands.cc with "EndOfEvent" argument<br>
 It have to pause the session command terminal.<br>
 Call SecondaryLoop to wait for exit event<br>
 @param aState
 @see : G4VisCommandReviewKeptEvents::SetNewValue
 */
void G4UIWt::PauseSessionStart (
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








void G4UIWt::ActivateCommand(
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
  printf("G4UIWt::ActivateCommand found : %s \n",targetCom.data());
#endif
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

void G4UIWt::InitHelpTreeAndVisParametersWidget()
{
  
  if (! fHelpTreeWidget ) {
    fHelpTreeWidget = new Wt::WTree();
  }
  
  // build widget
  fHelpTreeWidget->setSelectionMode(Wt::SingleSelection);
  
  fHelpTreeWidget->itemSelectionChanged ().connect(this,&G4UIWt::HelpTreeClicCallback);
}


/**   Fill the Help Tree Widget
 @param aParent : parent item to fill
 @param aCommandTree : commandTree node associate with this part of the Tree
 */
void G4UIWt::CreateHelpTree(
                                 Wt::WTreeNode *aParent
                                 ,G4UIcommandTree *aCommandTree
                                 )
{
  if (aParent == NULL) return;
  if (aCommandTree == NULL) return;
  
  
  // Creating new item
  Wt::WTreeNode * newItem;
  
  Wt::WString commandText = "";
  // Get the Sub directories
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {
    
    // trim path
    std::string whiteSpaces( " \f\n\r\t\v" );
    std::string path = (aCommandTree->GetTree(a+1)->GetPathName()).data();
    
    std::string::size_type posR = path.find_last_not_of( whiteSpaces );
    path.erase( posR + 1 );
    
    std::string::size_type posL = path.find_first_not_of( whiteSpaces );
    path.erase( 0, posL );
    
    commandText = Wt::WString(path.c_str());
    
    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,path.c_str());
    if (newItem == NULL) {
      newItem = new Wt::WTreeNode(GetShortCommandPath(path), 0,aParent);
    }
    CreateHelpTree(newItem,aCommandTree->GetTree(a+1));
  }
  
  // Get the Commands
  
  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    // trim path
    std::string whiteSpaces( " \f\n\r\t\v" );
    std::string path = (aCommandTree->GetCommand(a+1)->GetCommandPath()).data();
    
    std::string::size_type posR = path.find_last_not_of( whiteSpaces );
    path.erase( posR + 1 );
    
    std::string::size_type posL = path.find_first_not_of( whiteSpaces );
    path.erase( 0, posL );
    
    commandText = Wt::WString(path.c_str());
    
    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,path.c_str());
    if (newItem == NULL) {

      newItem = new Wt::WTreeNode(GetShortCommandPath(path), 0,aParent);
      
      newItem->collapse();
    }
  }
}







/** Find a treeItemWidget in the help tree
 @param aCommand item's String to look for
 @return item if found, NULL if not
 */
Wt::WTreeNode* G4UIWt::FindTreeItem(
                                                       Wt::WTreeNode *aParent
                                                       ,const std::string& aCommand
                                                       )
{
  if (aParent == NULL) return NULL;
  
  // Suppress last "/"
  std::string myCommand = aCommand;
  
  
  if (myCommand.rfind("/") == (myCommand.size()-1)) {
    myCommand = myCommand.substr(0,myCommand.size()-1);
  }
  
  if (GetLongCommandPath(aParent) == Wt::WString(myCommand.c_str()))
    return aParent;
  
  Wt::WTreeNode * tmp = NULL;
  for (unsigned int a=0;a<aParent->childNodes().size();a++) {
    if (!tmp)
      tmp = FindTreeItem(aParent->childNodes().at(a),myCommand);
  }
  return tmp;
}



/**   Build the command list parameters in a Wt::WString<br>
 Reimplement partialy the G4UIparameter.cc
 @param aCommand : command to list parameters
 @see G4UIparameter::List()
 @see G4UIcommand::List()
 @return the command list parameters, or "" if nothing
 */
Wt::WString G4UIWt::GetCommandList (
                                                 const G4UIcommand *aCommand
                                                 )
{
  
  Wt::WString txt ="";
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
    txt += "Command " + Wt::WString((char*)(commandPath).data()) + "\n";
  }
  txt += "Guidance :\n";
  
  for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
    txt += Wt::WString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
  }
  if( ! rangeString.isNull() ) {
    txt += " Range of parameters : " + Wt::WString((char*)(rangeString).data()) + "\n";
  }
  if( n_parameterEntry > 0 ) {
    G4UIparameter *param;
    
    // Re-implementation of G4UIparameter.cc
    
    for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ ) {
      param = aCommand->GetParameter(i_thParameter);
      txt += "\nParameter : " + Wt::WString((char*)(param->GetParameterName()).data()) + "\n";
      if( ! param->GetParameterGuidance().isNull() )
        txt += Wt::WString((char*)(param->GetParameterGuidance()).data())+ "\n" ;
      char myChar = param->GetParameterType();
      txt += " Parameter type  : " + Wt::WString(&myChar) + "\n";
      if(param->IsOmittable()){
        txt += " Omittable       : True\n";
      } else {
        txt += " Omittable       : False\n";
      }
      if( param->GetCurrentAsDefault() ) {
        txt += " Default value   : taken from the current value\n";
      } else if( ! param->GetDefaultValue().isNull() ) {
        txt += " Default value   : " + Wt::WString((char*)(param->GetDefaultValue()).data())+ "\n";
      }
      if( ! param->GetParameterRange().isNull() ) {
        txt += " Parameter range : " + Wt::WString((char*)(param->GetParameterRange()).data())+ "\n";
      }
      if( ! param->GetParameterCandidates().isNull() ) {
        txt += " Candidates      : " + Wt::WString((char*)(param->GetParameterCandidates()).data())+ "\n";
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
G4bool G4UIWt::IsGUICommand(
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
      char myChar = param->GetParameterType();
      if (myChar == 'd') {
        return true;
      }
      if (myChar == 'b') {
        return true;
      }
      if (myChar == 'i') {
        return true;
      }
      if (myChar == 's' && (!param->GetParameterCandidates().isNull())) {
        return true;
      }
    }
  }
  return false;
}


/**  Implement G4VBasicShell vurtual function
 */

G4bool G4UIWt::GetHelpChoice(
                                  G4int&
                                  )
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::GetHelpChoice SHOULD NEVER GO HERE");
#endif
  return true;
}





/***************************************************************************/
//
//             SLOTS DEFINITIONS
//
/***************************************************************************/

/**   Called when user give "help" command.
 */
void G4UIWt::ShowHelpCallback (
)
{
  TerminalHelp("");
}


/**   Called when user click on clear button. Clear the text Output area
 */
void G4UIWt::ClearButtonCallback (
)
{
  fCoutTBTextArea->setText("");
  fG4cout.removeRows(0,fG4cout.rowCount());
}

/**   Called when user exit session
 */
void G4UIWt::ExitSession (
)
{
  SessionTerminate();
}

void G4UIWt::ExitHelp(
) const
{
}




/**   Callback call when "enter" clicked on the command zone.<br>
 If command has no parameters :send the command to geant4
 Else, open a dialog for parameters input
 @param aCommand
 */
void G4UIWt::ButtonCallback (
                                  const char* aCommand
                                  )
{
  G4String ss = G4String(aCommand);
  ss = ss.strip(G4String::leading);
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();
  
  G4UIcommand* command = treeTop->FindPath(ss);
  
  if (command) {
    // if is GUI, then open a dialog
    if (IsGUICommand(command)) {
      Wt::WDialog* menuParameterDialog = new Wt::WDialog();
      
      if (CreateVisCommandGroupAndToolBox(command,menuParameterDialog,1,true)) {
        menuParameterDialog->setWindowTitle (aCommand);
        
        // exec this dialog, apply the command automaticaly, and return
        menuParameterDialog->exec();
        return;
      }
    }
  }
  
  ApplyShellCommand(ss,fExitSession,fExitPause);
  
  // Rebuild help tree
  FillHelpTree();
  
  if(fExitSession==true)
    SessionTerminate();
}





/**   Callback called when user give a new string to look for<br>
 Display a list of matching commands descriptions. If no string is set,
 will display the complete help tree
 */
void G4UIWt::LookForHelpStringCallback(
)
{
  Wt::WString searchText = fHelpLine->text();
  
  fHelpArea->setText("");
  if (searchText =="") {
    // clear old help tree
    fHelpTreeWidget = new Wt::WTree();
    
    FillHelpTree();
    
    return;
  } else {
    OpenHelpTreeOnCommand(searchText);
  }
}




Wt::WString G4UIWt::GetShortCommandPath(
                                                     const std::string& aTxt
                                                     )
{
  std::string commandPath;
  if (aTxt.find("/", 0) != std::string::npos) {
    //   commandPath = commandPath.right(commandPath.size()-1);
    commandPath = aTxt.substr(aTxt.size()-1);
  }
  
  //  commandPath = commandPath.right(commandPath.size()-commandPath.lastIndexOf("/",-2)-1);
  commandPath = commandPath.substr(commandPath.size()-commandPath.rfind("/",-2)-1);
  
  //  if (commandPath.lastIndexOf("/") == (commandPath.size()-1)) {
  if (commandPath.rfind("/") == (commandPath.size()-1)) {
    //    commandPath = commandPath.left(commandPath.size()-1);
    commandPath = commandPath.substr(0,commandPath.size()-1);
  }
  
  return commandPath.c_str();
}


Wt::WString G4UIWt::GetLongCommandPath(
                                                    Wt::WTreeNode* item
                                                    )
{
  if (item == NULL) return "";
  
  // rebuild path:
  Wt::WString itemText = "";
  itemText = item->label()->text();
  
  while (item->parentNode() != NULL) {
    itemText = item->parentNode()->label()->text()+"/"+itemText;
    item = item->parentNode();
  }
  itemText = "/"+itemText;
  
  return itemText;
}












G4WTabWidget::G4WTabWidget(
 Wt::WContainerWidget*& /* split */
):
 Wt::WTabWidget()
,tabSelected(false)
,lastCreated(-1)
{
}

G4WTabWidget::G4WTabWidget(
):Wt::WTabWidget()
,tabSelected(false)
,lastCreated(-1)
{
}



void G4UIWt::TabCloseCallback(int a){
  Wt::WWidget* temp =  fViewerTabWidget->widget(a);
  fViewerTabWidget->removeTab (temp);
  
  delete temp;
  
  if (fViewerTabWidget->count() == 0) {
    fEmptyViewerTabLabel->show();
  }
  
}


void G4UIWt::ToolBoxActivated(int a){
    
  if (fUITabWidget->widget(a) == fHelpTBWidget) {
    // Rebuild the help tree
    FillHelpTree();
  } else if (fUITabWidget->widget(a) == fSceneTreeComponentsTBWidget) {
    fSceneTreeComponentsTBWidget->hide();
  }
}


#endif
