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
// GEANT4 tag $Name: not supported by cvs2svn $
//
// L. Garnier

#ifdef G4UI_BUILD_WT_SESSION

#include "G4Types.hh"

#include <string.h>

#include "G4UIWt.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommandStatus.hh"

#include "G4Wt.hh"

#include <Wt/Wapplication>
#include <Wt/Wlineedit>
#include <Wt/Wwidget>
#include <Wt/Wlayout>
#include <Wt/Wpushbutton>
#include <Wt/Wlabel>
#include <Wt/Wscrollbar>
#include <Wt/Wdialog>
#include <Wt/Wevent>
#include <Wt/Wtextedit>
#include <Wt/Wtabwidget>
#include <Wt/WCanvasPaintDevice>
#include <Wt/WServer>

#include <Wt/WPaintedWidget>
#include <Wt/WPainter>

#include <Wt/WApplication>
#include <Wt/WVBoxLayout>
#include <Wt/WHBoxLayout>
#include <Wt/WPanel>
//#include <Wt/WSignalMapper>
#include <Wt/WText>
#include <Wt/Ext/MenuItem>
#include <Wt/WEnvironment>
#include <Wt/WGLWidget>
#include <Wt/WImage>

#include <Wt/Ext/Menu>
#include <Qlistwidget.h>


#include <stdlib.h>

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
// FOR TEST
class MyPaintedWidget : public Wt::WPaintedWidget
 {
 public:
   MyPaintedWidget(Wt::WContainerWidget *parent = 0)
     : Wt::WPaintedWidget(parent),
       foo_(100)
   {
      resize(200, 200); // provide a default size
   }

   void setFoo(int foo) {
      foo_ = foo;
      update(); // trigger a repaint
   }

 protected:
   void paintEvent(Wt::WPaintDevice *paintDevice) {
     Wt::WPainter painter(paintDevice);
     painter.drawLine(20, 20, foo_, foo_);
   }

 private:
   int foo_;
 };
// END FOR TEST





G4UIWt::G4UIWt (
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
,fLockTabUpdate(false)
{

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt creation 1\n");
#endif
  G4Wt* interactorManager = G4Wt::getInstance (argc,argv,(char*)"Wt");
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt creation 2\n");
#endif

  if (!interactorManager->GetMainInteractor()) {
    G4cout        << "G4UIWt : Unable to init Wt. Aborted" << G4endl;
  }
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt creation 3\n");
#endif
  
  UI = G4UImanager::GetUIpointer();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt creation 4\n");
#endif
  if(UI!=NULL) UI->SetSession(this);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt creation 5\n");
#endif
  if(UI!=NULL) UI->SetG4UIWindow(this);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt creation 6\n");
#endif

  CreateWidgets ();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt END :%d\n",this);
#endif
}


void G4UIWt::CreateWidgets ()
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets\n");
#endif
  // Check if already define in external app WMainWindow
  bool found = false;
//   foreach (WWidget *widget, WApplication::allWidgets()) {
//     if ((found== false) && (widget->inherits("WMainWindow"))) {
//       found = true;
//     }
//   }

  if (found) {
    G4cout        << "G4UIWt : Found an external App with a WMainWindow already defined. Aborted" << G4endl;
    return;
  }

  fMainWindow = new Wt::WContainerWidget(wApp->root());

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::Initialise after main window creation +++++++++++:%d\n",wApp->root());
#endif

  //  Wt::WWidget *mainWidget = new Wt::WWidget(fMainWindow);

  Wt::WContainerWidget *mainWidget = new Wt::WContainerWidget(fMainWindow);
  Wt::WVBoxLayout* mainWindowVLayout = new Wt::WVBoxLayout();
  fMainWindow->setLayout(mainWindowVLayout);
  fMainWindow->resize(1400,800);

  fMyVSplitter = new Wt::Ext::Splitter(Wt::Horizontal,fMainWindow);
  fMyVSplitter->setOrientation(Wt::Horizontal);
  fToolBox = new Wt::WContainerWidget();

  // Set layouts


  Wt::WContainerWidget* commandLineWidget = new Wt::WContainerWidget();
  Wt::WVBoxLayout *layoutCommandLine = new Wt::WVBoxLayout();


#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt :: 5\n");
#endif

  //  mainWidget->setLayout(mainWindowVLayout);

  //  fMainWindow->setCentralWidget(mainWidget);


#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt :: 6\n");
#endif
  // Add a quit subMenu
  Wt::WContainerWidget *ex = new Wt::WContainerWidget();
  
  fToolBar = new Wt::Ext::ToolBar(ex);
  Wt::Ext::Menu *fileMenu = new Wt::Ext::Menu();
  fileMenu->addItem("File",this, &G4UIWt::ExitSession);

  fToolBar->addButton("File", fileMenu);
  mainWindowVLayout->addWidget(ex);




  AddInteractor ("file",(G4Interactor)fileMenu);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt :: 7\n");
#endif


  // fill them

  fCommandLabel = new Wt::WLabel("Command:");

  fCommandArea = new Wt::WLineEdit();
  fCommandArea->keyPressed().connect(this,&G4UIWt::CommandLineSlot);
  //  fCommandArea->installEventFilter(this);
  //  fCommandArea->activateWindow();




  //  fCommandArea->setFocusPolicy ( Wt::StrongFocus );
  //  fCommandArea->setFocus(Wt::TabFocusReason);

  layoutCommandLine->addWidget(fCommandLabel);
  layoutCommandLine->addWidget(fCommandArea);


  fHelpTBWidget = new Wt::WPanel();
  //  fHelpTBWidget->setCollapsed(true);
  fHelpTBWidget->setCollapsible(true);
  fHistoryTBWidget = new Wt::WPanel();
  fHistoryTBWidget->setCollapsed(true);
  fHistoryTBWidget->setCollapsible(true);
  fCoutTBWidget = new Wt::WPanel();
  fCoutTBWidget->setCollapsible(true);
  fCoutTBWidget->setCollapsed(true);

  fVisParametersTBWidget = new Wt::WPanel();
  fVisParametersTBWidget->setCollapsed(true);
  fVisParametersTBWidget->setCollapsible(true);
  fViewComponentsTBWidget = new Wt::WPanel();
  fViewComponentsTBWidget->setCollapsed(true);
  fViewComponentsTBWidget->setCollapsible(true);
  
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets 7\n");
#endif
  CreateVisParametersTBWidget();

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets 8\n");
#endif
  CreateViewComponentsTBWidget();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets 9\n");
#endif
  CreateHelpTBWidget();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets 10\n");
#endif
  CreateCoutTBWidget();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets 11\n");
#endif
  CreateHistoryTBWidget();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets 12\n");
#endif

  // the splitter 
  //  fToolBox->addItem(fVisParametersTBWidget,"Vis parameters");
  //  fToolBox->addItem(fViewComponentsTBWidget,"Viewer components");
  fToolBox->addWidget(fHelpTBWidget);
  fToolBox->addWidget(fCoutTBWidget);
  fToolBox->addWidget(fHistoryTBWidget);

  fHelpTBWidget->setTitle("Help");
  fCoutTBWidget->setTitle("Cout");
  fHistoryTBWidget->setTitle("History");


  //  fToolBox->setSizePolicy (WSizePolicy(WSizePolicy::Fixed,WSizePolicy::Fixed));

  fEmptyViewerTabLabel = new Wt::WLabel("         If you want to have a Viewer, please use /vis/open commands. ");

  // FIXME : For test
#define TEST
#ifdef TEST
   fTabWidget = new G4WTabWidget();
#else
  G4WTabWidget *fTabWidget = new G4WTabWidget();
#endif
  fTabWidget->resize(300,300);
  fEmptyViewerTabLabel->show();
  fTabWidget->show();

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CreateWidgets   tabWidget as Id:%d\n",fTabWidget);
#endif
  Wt::WContainerWidget *myTabCont = new Wt::WContainerWidget();
  
  //  fTabWidget->addWidget(myTab);

  // GL
  myTabCont->resize(220,220);
  fEmptyViewerTabLabel->resize(200,200);

  myTabCont->addWidget(fEmptyViewerTabLabel);
  fTabWidget->addTab(myTabCont,"titi");
  MyPaintedWidget *e = new MyPaintedWidget(myTabCont);
  // FIXME : For test


  // Only at creation. Will be set visible when sessionStart();
  fEmptyViewerTabLabel->hide();


  fToolBox->resize(400,800);
  fEmptyViewerTabLabel->resize(400,400);
  fMyVSplitter->addWidget(fToolBox);
  fMyVSplitter->addWidget(fEmptyViewerTabLabel);
  fMyVSplitter->addWidget(fTabWidget);

  commandLineWidget->setLayout(layoutCommandLine);
  //  commandLineWidget->setSizePolicy (WSizePolicy(WSizePolicy::Minimum,WSizePolicy::Minimum));
  mainWindowVLayout->addWidget(fMyVSplitter,1);
  mainWindowVLayout->addWidget(commandLineWidget);




  // Connect signal
  fCommandArea->enterPressed().connect(this,&G4UIWt::CommandEnteredCallback);
  fHelpTBWidget->expanded().connect(this,&G4UIWt::HelpToolBoxActivated);

  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP


  //  fMainWindow->setWindowTitle( tr("G4UI Session") ); 
  //  fMainWindow->resize(900,600); 
  //  fMainWindow->move(WPoint(50,100));

  // Set not visible until session start
  fMainWindow->hide();

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt END\n");
#endif
  // FIXME::
  //  AddTabWidget(new Wt::WLabel("Premier..."),"my name",50,50);
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
void G4UIWt::CreateHistoryTBWidget(
) 
{

  //  Wt::WVBoxLayout *layoutHistoryTB = new Wt::WVBoxLayout();
  fHistoryTBTableList = new Wt::WSelectionBox();
  fHistoryTBTableList->setSelectionMode(Wt::SingleSelection);
  fHistoryTBTableList->changed().connect(this,&G4UIWt::CommandHistoryCallback);
  fHistoryTBTableList->keyPressed ().connect(this,&G4UIWt::ChangeFocus);
  fHistoryTBTableList->resize(350,250);
  //  fHistoryTBTableList->installEventFilter(this);

  //  layoutHistoryTB->addWidget(fHistoryTBTableList);

  fHistoryTBWidget->setCentralWidget(fHistoryTBTableList);
}

/** Create the Help ToolBox Widget
 */
void G4UIWt::CreateHelpTBWidget(
) 
{

  
  Wt::WContainerWidget *helpWidget = new Wt::WContainerWidget();
  Wt::WHBoxLayout *helpLayout = new Wt::WHBoxLayout();
  Wt::WVBoxLayout *vLayout = new Wt::WVBoxLayout();
  fHelpVSplitter = new Wt::Ext::Splitter(Wt::Horizontal);


  fHelpLine = new Wt::WLineEdit();
  fHelpTBWidget->setCentralWidget(helpWidget);	
  helpWidget->addWidget(new Wt::WLabel(Wt::WString("Search :")));
  helpWidget->addWidget(fHelpLine);
  fHelpLine->enterPressed().connect(this, &G4UIWt::LookForHelpStringCallback );
  
  // Create Help tree
  FillHelpTree();
  return;
  
  fHelpArea = new Wt::WTextEdit(fHelpVSplitter);
  fHelpArea->setReadOnly(true);
  
  // Set layouts
  
  if (fHelpTreeWidget) {
    fHelpVSplitter->addWidget(fHelpTreeWidget);
  }
  fHelpVSplitter->addWidget(fHelpArea);
  
  
  vLayout->addWidget(helpWidget);
  vLayout->addWidget(fHelpVSplitter,1);
  
  fHelpTBWidget->setMinimumSize(50,50);
  //  fHelpTBWidget->setSizePolicy (WSizePolicy(WSizePolicy::Minimum,WSizePolicy::Minimum));
  // set the splitter size
  //  Wt::WList<int> list;
//   list.append( 50 );
//   list.append( 50 );
//  fHelpVSplitter->setSizes(list);
  
  helpWidget->setLayout(helpLayout);
  fHelpTBWidget->setCentralWidget(helpWidget);
}


/** Create the Cout ToolBox Widget
 */
void G4UIWt::CreateCoutTBWidget(
) 
{


  Wt::WContainerWidget* panelContainer = new Wt::WContainerWidget();
  Wt::WVBoxLayout *panelContainerVLayout = new Wt::WVBoxLayout();
  panelContainer->setLayout(panelContainerVLayout);

  fCoutTBTextArea = new Wt::WTextArea();
  fCoutFilter = new Wt::WLineEdit();
  Wt::WLabel* coutFilterLabel = new Wt::WLabel("Filter : ");
  panelContainerVLayout->addWidget(fCoutTBTextArea);
  panelContainerVLayout->addWidget(fCoutFilter);
  panelContainerVLayout->addWidget(coutFilterLabel);
  
  Wt::WPushButton *coutTBClearButton = new Wt::WPushButton("clear",panelContainer);
  coutTBClearButton->clicked().connect(this,&G4UIWt::ClearButtonCallback);

  //  fCoutFilter->changed().connect(this, (&G4UIWt::CoutFilterCallback));

  Wt::WSignalMapper<Wt::WString> *mySignalMapper = new Wt::WSignalMapper<Wt::WString>(this);

  //MenuItem::activated()

  // Following line cause :
  // _____________________
  // dyld: lazy symbol binding failed: Symbol not found: __ZN5boost7signals6detail11signal_baseC2ERKNS_9function2IbNS1_12stored_groupES4_SaINS_13function_baseEEEERKNS_3anyE
  //   Referenced from: /Users/garnier/Work/geant4/lib/Darwin-g++/libG4UIbasic.dylib
  //   Expected in: flat namespace
  
  // dyld: Symbol not found: __ZN5boost7signals6detail11signal_baseC2ERKNS_9function2IbNS1_12stored_groupES4_SaINS_13function_baseEEEERKNS_3anyE
  //   Referenced from: /Users/garnier/Work/geant4/lib/Darwin-g++/libG4UIbasic.dylib
  //   Expected in: flat namespace


  // FIXME WT
  mySignalMapper->mapped().connect(this, &G4UIWt::CoutFilterCallback);

  fCoutFilter->changed().connect(mySignalMapper, &Wt::WSignalMapper<Wt::WString>::map);
  mySignalMapper->setMapping(fCoutFilter, fCoutFilter->text());
  


  fCoutTBTextArea->setReadOnly(true);

  Wt::WContainerWidget* coutButtonWidget = new Wt::WContainerWidget(panelContainer);
  Wt::WHBoxLayout* layoutCoutTBButtons = new Wt::WHBoxLayout(coutButtonWidget);
  layoutCoutTBButtons->addWidget(coutTBClearButton);
  layoutCoutTBButtons->addWidget(coutFilterLabel);
  layoutCoutTBButtons->addWidget(fCoutFilter);

  panelContainerVLayout->addWidget(coutButtonWidget);

  fCoutTBWidget->setCentralWidget(panelContainer);
}


/** Create the VisParameters ToolBox Widget
 */
void G4UIWt::CreateVisParametersTBWidget(
) 
{
}


/** Create the ViewComponents ToolBox Widget
 */
void G4UIWt::CreateViewComponentsTBWidget(
) 
{
}


/**   Add a new tab widget.
  Create the tab if it was not done
*/
bool G4UIWt::AddTabWidget(
 Wt::WWidget* aWidget
 ,Wt::WString name
,int sizeX
,int sizeY
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::AddTabWidget %d %d\n",sizeX, sizeY);
#endif

  if (fTabWidget == NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget +++++\n");
#endif
    fTabWidget = new G4WTabWidget();
    fTabWidget->resize(200,200);
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget 1\n");
#endif
    //    fTabWidget->setTabCloseable (true); 
    
    //    fTabWidget->setUsesScrollButtons (true);
    
    //    fTabWidget->setSizePolicy (WSizePolicy(WSizePolicy::Maximum,WSizePolicy::Maximum));
    
    //    Wt::WSizePolicy policy = fTabWidget->sizePolicy();
    //    policy.setHorizontalStretch(1);
    //    policy.setVerticalStretch(1);
    //    fTabWidget->setSizePolicy(policy);
    
    fTabWidget->tabClosed().connect(this, &G4UIWt::TabCloseCallback);
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget 2\n");
#endif
    fTabWidget->currentChanged().connect(this,&G4UIWt::UpdateTabWidget); 
  }

  fLastWTabSizeX = sizeX;
  fLastWTabSizeY = sizeY;

  if (!aWidget) {
    return false;
  }

  // Remove WLabel 

  // L.Garnier 26/05/2010 : not exactly the same in qt3. Could cause some
  // troubles
  if (fEmptyViewerTabLabel != NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget 4\n");
#endif
    if ( fMyVSplitter->isVisible() != -1) {
      
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::AddTabWidget 5\n");
#endif
      fMyVSplitter->removeWidget (fEmptyViewerTabLabel);
      //    fEmptyViewerTabLabel->setParent(NULL);
      delete fEmptyViewerTabLabel;
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::AddTabWidget 5a\n");
#endif
      fEmptyViewerTabLabel = NULL;
      
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::AddTabWidget 5a2 %d %d\n",fMyVSplitter,fTabWidget);
#endif
      //      fMyVSplitter->addWidget(fTabWidget);
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::AddTabWidget 5b\n");
#endif
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::AddTabWidget 5c\n");
#endif
      fTabWidget->show();
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::AddTabWidget 6\n");
#endif
      
      // FIXME :
      // Not necessary ?
      //    aWidget->setParent(fTabWidget);
    }
  }
  
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::AddTabWidget ADD %d %d + %d %d---------------------------------------------------\n",sizeX, sizeY,sizeX-fTabWidget->width().value(),sizeY-fTabWidget->height().value());
#endif

  if (fMainWindow->isVisible()) {

    // get the size of the tabbar
    double tabBarX = 0;
    double tabBarY = 0;
    if (fTabWidget->count() >0) {
      tabBarX = fTabWidget->width().value()-fTabWidget->widget(0)->width().value();
      tabBarY = fTabWidget->height().value()-fTabWidget->widget(0)->height().value();
    }

    //    fMainWindow->resize(tabBarX+fMainWindow->width().value()+sizeX-fTabWidget->width().value(),tabBarY+fMainWindow->height().value()+sizeY-fTabWidget->height().value());
  }

  // Problems with resize. The widgets are not realy drawn at this step,
  // then we have to force them on order to check the size

#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget 7\n");
#endif
  fLastTabWidgetInsert = new Wt::WContainerWidget();
  fTabWidget->addTab(fLastTabWidgetInsert,name);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::AddTabWidget 8 tabWidget is now :%d   containerId is :%d\n",fTabWidget,fLastTabWidgetInsert);
#endif
  
  fTabWidget->setTabCloseable (fTabWidget->currentIndex(),true); 

  fTabWidget->setCurrentIndex(fTabWidget->count()-1);

#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget 9\n");
#endif

  // Set visible
   fTabWidget->setCurrentIndex(fTabWidget->count()-1);
  
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::AddTabWidget END\n");
#endif
  return true;
}


void G4UIWt::UpdateTabWidget(int tabNumber) {
  if (fLockTabUpdate)
    return;
  fLockTabUpdate = true;

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::UpdateTabWidget %d\n",tabNumber);
#endif
  if (  fTabWidget == NULL) {
    fTabWidget = new G4WTabWidget;
  }
  

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::UpdateTabWidget CALL REPAINT tabGL\n");
#endif

  fTabWidget->setCurrentIndex(tabNumber);

  // Send this signal to unblock graphic updates !
  fTabWidget->setTabSelected(false);

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::UpdateTabWidget 1\n");
#endif
  //  fTabWidget->show();

  // This will send a paintEvent to OGL Viewers
  fTabWidget->setTabSelected(true);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::UpdateTabWidget 2\n");
#endif

  // FIXME :
  //  Wt::WCoreApplication::sendPostedEvents () ;

  fLockTabUpdate = false;
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::UpdateTabWidget END\n");
#endif
}


/** Send resize event to all tabs
 */
// void G4UIWt::ResizeTabWidget( WResizeEvent* e) {
//   for (G4int a=0;a<fTabWidget->count() ;a++) {
// #ifdef G4DEBUG_INTERFACES_BASIC
//     printf("G4UIWt::ResizeTabWidget +++++++++++++++++++++++++++++++++++++++\n");
// #endif
//     fTabWidget->widget(a)->resize(e->size());
//   }
// }


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

  if (fEmptyViewerTabLabel != NULL) {
    bool visible = false;
    if (fTabWidget != NULL) {
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4UIWt::G4UIWt SessionStart 1a:%d\n",fMainWindow);
#endif
      fEmptyViewerTabLabel->show();
      fTabWidget->show();
    } else {
      fEmptyViewerTabLabel->show();
    }
  }

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 2:%d\n",fMainWindow);
#endif
  fMainWindow->show();
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 2b\n");
#endif
  // get the size of the tabbar
  double tabBarX = 0;
  double tabBarY = 0;

  if (fTabWidget != NULL) {
    tabBarX = -fTabWidget->widget(0)->width().value();
    tabBarY = -fTabWidget->widget(0)->height().value();
  }
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 2c\n");
#endif
  fMainWindow->resize(tabBarX+fMainWindow->width().value()+fLastWTabSizeX,tabBarY+fMainWindow->height().value()+fLastWTabSizeY);

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 3\n");
#endif
  // FIXME :
  //  Wt::WCoreApplication::sendPostedEvents () ;

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 4\n");
#endif
  //  interactorManager->DisableSecondaryLoop (); // TO KEEP

  // add a single entry point, at the default location (as determined
  // by the server configuration's deploy-path)

 
  //  if ((WApplication*)interactorManager->GetMainInteractor())
  //    ((WApplication*)interactorManager->GetMainInteractor())->exec();

  // on ne passe pas le dessous ? FIXME ????
  // je ne pense pas 13/06

  //   void* event; // TO KEEP
  //   while((event = interactorManager->GetEvent())!=NULL) {  // TO KEEP
  //     interactorManager->DispatchEvent(event); // TO KEEP
  //     if(exitSession==true) break; // TO KEEP
  //   } // TO KEEP

  /*
    interactorManager->EnableSecondaryLoop ();
  */
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::G4UIWt SessionStart 5\n");
#endif
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

  fCommandLabel->setText((char*)aPrompt.data());
}



void G4UIWt::SessionTerminate (
)
{
  G4Wt* interactorManager = G4Wt::getInstance ();
  //  quit();
  wApp->quit(); 
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
void G4UIWt::SecondaryLoop (
 G4String aPrompt
)
{
  if (!aPrompt) return;

  G4Wt* interactorManager = G4Wt::getInstance (); // TO KEEP ?
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
G4int G4UIWt::ReceiveG4cout (
 G4String aString
 )
{
  if (!aString) return 0;
  
  QStringList newStr;
  
  // Add to stringList
  newStr = QStringList(QString((char*)aString.data()).trimmed());
  fG4cout += newStr;
 
  QStringList result = newStr.filter(QString((fCoutFilter->text().toUTF8().c_str())));

  if (result.join("\n").isEmpty()) {
    return 0;
  }
  fCoutTBTextArea->setText(fCoutTBTextArea->text()+"\n"+result.join("\n").toStdString());
  fCoutTBTextArea->refresh();

  //  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());

  return 0;
}


/**
   Receive a cerr from Geant4. We have to display it in the cout zone
   @param aString : label to add in the display area
   @return 0
*/
G4int G4UIWt::ReceiveG4cerr (
 G4String aString
)
{
  if (!aString) return 0;

  QStringList newStr;

  // Add to stringList
  newStr = QStringList(QString((char*)aString.data()).trimmed());
  fG4cout += newStr;
 
  QStringList result = newStr.filter(QString((fCoutFilter->text().toUTF8().c_str())));

  fCoutTBTextArea->setText(fCoutTBTextArea->text()+"<font color=red>"+result.join("\n").toStdString()+"</font>");
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::ReceiveG4cerr : %s\n",result.join("\n").toStdString().c_str());
#endif
  //  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());
  fCoutTBTextArea->refresh();
  return 0;
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

  Wt::Ext::Menu *fileMenu = new Wt::Ext::Menu();
  fToolBar->addButton(aLabel,fileMenu); 

  AddInteractor (aName,(G4Interactor)fileMenu);
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

  Wt::Ext::Menu *parent = (Wt::Ext::Menu*)GetInteractor(aMenu);

  if(parent==NULL) return;
  
  Wt::WSignalMapper<Wt::WString> *signalMapper = new Wt::WSignalMapper<Wt::WString>(this);
  //MenuItem::activated()
  Wt::Ext::MenuItem *myMenu = new Wt::Ext::MenuItem(aLabel);
  signalMapper->mapped().connect(this, &G4UIWt::ButtonCallback);

  myMenu->activated().connect(signalMapper, &Wt::WSignalMapper<Wt::WString>::map);
  signalMapper->setMapping(myMenu, Wt::WString(aCommand));

  parent->add(myMenu);
//  signalMapper->mapConnect(myMenu->activated(), Wt::WString(aCommand));
//Wt::WAction *action = parent->addAction(aLabel, signalMapper, &Wt::WSignalMapper::map());

//  signalMapper->setMapping(action, );
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

  fHelpTBWidget->setCollapsed(false);
  fHistoryTBWidget->setCollapsed(true);
  fVisParametersTBWidget->setCollapsed(true);
  fViewComponentsTBWidget->setCollapsed(true);
  fCoutTBWidget->setCollapsed(true);
}



/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIWt::InitHelpTree()
{

#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 1\n");
#endif
  if (! fHelpTreeWidget ) {
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 1b\n");
#endif
    fHelpTreeWidget = new Wt::WTreeTable ();
    fHelpTreeRoot = new Wt::WTree();
  }
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 2\n");
#endif


  // build widget
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 3\n");
#endif
  fHelpTreeRoot->setSelectionMode(Wt::SingleSelection);
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 3b\n");
#endif
  return;
  //  fHelpTreeWidget->addColumn("Command");
  fHelpTreeWidget->setTree(fHelpTreeRoot, "Command");
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 4\n");
#endif



#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() 5\n");
#endif
  fHelpTreeRoot->itemSelectionChanged ().connect(this,&G4UIWt::HelpTreeClicCallback);  
  //  connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked (WTreeNode*,int)),this, &G4UIWt::HelpTreeDoubleClicCallback());  
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::InitHelpTree() END\n");
#endif

}
/**
   Create the help tree widget
   @param parent : parent of tree widget
   @return the widget containing the tree or NULL if it could not have beeen created
 */

void G4UIWt::FillHelpTree()
{
  //FIXME 
  if (! fHelpTreeWidget ) {
    InitHelpTree();
  }
  return; 

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

    //FIXME
    commandText = Wt::WString((char*)(treeTop->GetTree(a+1)->GetPathName()).data());

    // if already exist, don't create it !
    for (unsigned int b=0;b<fHelpTreeRoot->treeRoot()->childNodes().size() ;b++) {
      if (!newItem)
        newItem = FindTreeItem(fHelpTreeRoot->treeRoot()->childNodes()[b],commandText);
    }

    if (newItem == NULL) {
      
      newItem = new Wt::WTreeNode(GetShortCommandPath(commandText));
      fHelpTreeRoot->treeRoot()->addChildNode(newItem);
      //      newItem->setText(0,);
    }

    // look for childs
    CreateChildTree(newItem,treeTop->GetTree(a+1));
  }

}



/**   Fill the Help Tree Widget
   @param aParent : parent item to fill
   @param aCommandTree : commandTree node associate with this part of the Tree
*/
void G4UIWt::CreateChildTree(
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

    //    commandText = Wt::WString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).trimmed();
    commandText = Wt::WString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data());
    
    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
      newItem = new Wt::WTreeNode(GetShortCommandPath(commandText));
      aParent->addChildNode(newItem);
    }
    CreateChildTree(newItem,aCommandTree->GetTree(a+1));
  }



  // Get the Commands

  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    QStringList stringList;
    //    commandText = Wt::WString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed();
    commandText = Wt::WString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data());

    // if already exist, don't create it !
    newItem = FindTreeItem(aParent,commandText);
    if (newItem == NULL) {
      newItem = new Wt::WTreeNode(GetShortCommandPath(commandText));
      aParent->addChildNode(newItem);
      newItem->collapsed();
    }
  }
}

 
/** Find a treeItemWidget in the help tree
    @param aCommand item's String to look for
    @return item if found, NULL if not
*/
 Wt::WTreeNode* G4UIWt::FindTreeItem(
 Wt::WTreeNode *aParent
,const Wt::WString& aCommand
)
{
  if (aParent == NULL) return NULL;

  // Suppress last "/"
  QString myCommandQ = aCommand.toUTF8().c_str();
  
  if (myCommandQ.lastIndexOf("/") == (myCommandQ.size()-1)) {
    myCommandQ = myCommandQ.left(myCommandQ.size()-1);
  }
  Wt::WString myCommand = myCommandQ.toStdString().c_str();

  if (GetLongCommandPath(aParent) == myCommand)
    return aParent;
  
  Wt::WTreeNode * tmp = NULL;
  for (unsigned int a=0;a<aParent->childNodes().size();a++) {
    if (!tmp)
      tmp = FindTreeItem(aParent->childNodes()[a],myCommand);
  }
  return tmp;
}



/**   Build the command list parameters in a WString<br>
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
      txt += " Parameter type  : " + Wt::WString(QString(QChar(param->GetParameterType())).toStdString().c_str()) + "\n";
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


/**   Event filter method. Every event from WtApplication goes here.<br/>
   We apply a filter only for the Up and Down Arrow press when the WLineEdit<br/>
   is active. If this filter match, Up arrow we give the previous command<br/>
   and Down arrow will give the next if exist.<br/>
   @param obj Emitter of the event
   @param event Kind of event
*/
#ifdef QT
bool G4UIWt::eventFilter( // Should stay with a minuscule eventFilter because of Wt
 Wt::WObject *aObj
,WEvent *aEvent
)
{
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
      Wt::WKeyEvent *e = static_cast<WKeyEvent*>(aEvent);
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
        G4String ss = Complete(fCommandArea->text().toStdString().c_str());
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
  bool res= false;
  // change cursor position if needed
  if (moveCommandCursor == true) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::eventFilter setCursor Position\n");
#endif
    fCommandArea->setCursorPosition ( fCommandArea->text().length() );
    fCommandArea->setCursorPosition (4);
  } else {
    // pass the event on to the parent class
    res = Wt::WObject::eventFilter(aObj, aEvent);
  }
  return res;
}

#endif




/***************************************************************************/
//
//             SLOTS DEFINITIONS
//
/***************************************************************************/

void G4UIWt::CommandLineSlot(Wt::WKeyEvent e) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandLineSlot\n");
#endif

  bool moveCommandCursor = false;
  if ((e.key() == (Wt::Key_Down)) ||
      (e.key() == (Wt::Key_PageDown)) ||
      (e.key() == (Wt::Key_Up)) ||
      (e.key() == (Wt::Key_PageUp))) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandLineSlot Arrow\n");
#endif
    int selection = fHistoryTBTableList->currentIndex();
    if (fHistoryTBTableList->count()) {
      if (selection == -1) {
        selection = fHistoryTBTableList->count()-1;
      } else {
        if (e.key() == (Wt::Key_Down)) {
          if (selection <(fHistoryTBTableList->count()-1))
            selection++;
        } else if (e.key() == (Wt::Key_PageDown)) {
          selection = fHistoryTBTableList->count()-1;
        } else if (e.key() == (Wt::Key_Up)) {
          if (selection >0)
            selection --;
        } else if (e.key() == (Wt::Key_PageUp)) {
          selection = 0;
        }
      }
      //        fHistoryTBTableList->clearSelection();
      fHistoryTBTableList->setCurrentIndex(selection);
      //        fHistoryTBTableList->setCurrentItem(fHistoryTBTableList->item(selection));
    }
    moveCommandCursor = true;
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::eventFilter setCursor Position\n");
#endif
  } else if (e.key() == (Wt::Key_Tab)) {
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandLineSlot Tab\n");
#endif
    G4String ss = Complete(fCommandArea->text().toUTF8());
    fCommandArea->setText((char*)(ss.data()));

    // do not pass by parent, it will disable widget tab focus !
    
    // L.Garnier : MetaModifier is CTRL for MAC, but I don't want to put a MAC 
    // specific #ifdef
//   } else if (((e.modifiers () == Wt::ControlModifier) || (e.modifiers () == Wt::MetaModifier)) && (e.key() == Wt::Key_A)) {
//     fCommandArea->home(false);
//   } else if (((e.modifiers () == Wt::ControlModifier) || (e.modifiers () == Wt::MetaModifier)) && (e.key() == Wt::Key_E)) {
//       fCommandArea->end(false);
  }

}

void G4UIWt::ChangeFocus () {
  fCommandArea->setFocus();
}

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
  fG4cout.clear();
}

/**   Called when user exit session
*/
void G4UIWt::ExitSession (
)
{
  SessionTerminate();
}

void G4UIWt::ExitHelp(
)
{
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
    fHistoryTBTableList->addItem(fCommandArea->text());
    fHistoryTBTableList->clearSelection();
    fHistoryTBTableList->setCurrentIndex(0);
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


/**   Callback call when "enter" clicked on the command zone.<br>
   Send the command to geant4
   @param aCommand
*/
void G4UIWt::ButtonCallback (
 const Wt::WString& aCommand
)
{
  G4String ss = G4String(aCommand.toUTF8().c_str());
  ApplyShellCommand(ss,exitSession,exitPause);

  // Rebuild help tree
  FillHelpTree();

  if(exitSession==true) 
    SessionTerminate();
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
  
  std::set< Wt::WTreeNode * > list =fHelpTreeRoot->selectedNodes();
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
 
#ifdef QT
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
  
  Wt::WList<WTreeNode *> list =fHelpTreeWidget->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
  if (!item)
    return;

  fCommandArea->clear();
  fCommandArea->setText(GetLongCommandPath(item));
}
#endif


/**   Callback called when user select an old command in the command history<br>
   Give it to the command area.
*/
void G4UIWt::CommandHistoryCallback(
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
  printf("G4UIWt::CommandHistoryCallback\n");
#endif
  //  int item = -1;
  if (!fHistoryTBTableList)
    return ;

  
//   std::set< int >& list =fHistoryTBTableList->selectedIndexes ();
//   if (list.empty())
//     return;
//   item = *list.begin();
//   if (item < 0)
//     return;
  fCommandArea->setText(fHistoryTBTableList->currentText ());
#ifdef G4DEBUG_INTERFACES_BASIC
    printf("G4UIWt::CommandHistoryCallback change text\n");
#endif
}


void G4UIWt::CoutFilterCallback(
  //const Wt::WString & text) {
Wt::WString text) {

  QStringList result = fG4cout.filter(text.toUTF8().c_str());
  fCoutTBTextArea->setText(result.join("\n").toStdString());

  //  fCoutTBTextArea->repaint();
  //  fCoutTBTextArea->verticalScrollBar()->setSliderPosition(fCoutTBTextArea->verticalScrollBar()->maximum());

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
    fHelpTreeRoot = new Wt::WTree();

    FillHelpTree();

    return;
  } else {
    OpenHelpTreeOnCommand(searchText);
  }
}


void G4UIWt::OpenHelpTreeOnCommand(
 const Wt::WString & searchText
)
{
//_________________________________________________________________________________
#ifdef XXXXXXX

  // the help tree
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4UIcommandTree * treeTop = UI->GetTree();
  
  G4int treeSize = treeTop->GetTreeEntry();

  // clear old help tree
  fHelpTreeRoot = new Wt::WTree();

  // look for new items

  int tmp = 0;

  QMap<int,QString> commandResultMap;
  QMap<int,QString> commandChildResultMap;

  for (int a=0;a<treeSize;a++) {
    G4UIcommand* command = treeTop->FindPath(treeTop->GetTree(a+1)->GetPathName().data());
    tmp = QString(GetCommandList (command).toUTF8().c_str()).count(QString(searchText.toUTF8().c_str()),Qt::CaseInsensitive);
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
  fHelpTreeRoot->setSelectionMode(Wt::SingleSelection);
  fHelpTreeWidget->setTree(fHelpTreeRoot, Wt::WString("Command"));
  fHelpTreeWidget->addColumn(Wt::WString("Match"),0);
 
 //  fHelpTreeWidget->setColumnCount(2);
  //  QStringList labels;
  //  labels << Qt::WString("Command") << Wt::WString("Match");
  //  fHelpTreeWidget->setHeaderLabels(labels);

  if (commandResultMap.empty()) {
    fHelpArea->setText("No match found");
    return;
  }

  QMap<int,QString>::const_iterator i = commandResultMap.constEnd();
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
    //    newItem = new Wt::WTreeNode(fHelpTreeRoot);
    //    Wt::WString commandStr = i.value().trimmed();
    QString commandStr = i.value();

    if (commandStr.indexOf("/") == 0) {
      commandStr = commandStr.right(commandStr.size()-1);
    }
      
    newItem = new Wt::WTreeNode(commandStr.toStdString());
    //    newItem->setText(0,commandStr.toStdString());
    newItem->setText(1,progressStr);
    
    newItem->setForeground ( 1, Wt::WBrush(Wt::blue) );
    progressStr = "|";
    i--;
  }
  fHelpTreeWidget->resizeColumnToContents (0);
  fHelpTreeWidget->sortItems(1,Wt::DescendingOrder);
  //  fHelpTreeWidget->setColumnWidth(1,10);//resizeColumnToContents (1);


#endif
//_________________________________________________________________________________
}




QMap<int,QString> G4UIWt::LookForHelpStringInChildTree(
 G4UIcommandTree *aCommandTree
,const Wt::WString & text
 )
{
  QMap<int,QString> commandResultMap;
  if (aCommandTree == NULL) return commandResultMap;
  

  // Get the Sub directories
  int tmp = 0;
  QMap<int,QString> commandChildResultMap;
  
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {
    const G4UIcommand* command = aCommandTree->GetGuidance();
    tmp = QString(GetCommandList (command).toUTF8().c_str()).count(QString(text.toUTF8().c_str()),Qt::CaseInsensitive);
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
    tmp = QString(GetCommandList (command).toUTF8().c_str()).count(QString(text.toUTF8().c_str()),Qt::CaseInsensitive);
    if (tmp >0) {
      commandResultMap.insertMulti(tmp,QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()));
    }
    
  }
  return commandResultMap;
}

  
Wt::WString G4UIWt::GetShortCommandPath(
Wt::WString commandPath
)
{
  QString commandPathQ = commandPath.toUTF8().c_str();

  if (commandPathQ.indexOf("/") == 0) {
    commandPathQ = commandPathQ.right(commandPathQ.size()-1);
  }

  commandPathQ = commandPathQ.right(commandPathQ.size()-commandPathQ.lastIndexOf("/",-2)-1);
 
 if (commandPathQ.lastIndexOf("/") == (commandPathQ.size()-1)) {
    commandPathQ = commandPathQ.left(commandPathQ.size()-1);
 }

 return Wt::WString(commandPathQ.toStdString());
}


Wt::WString G4UIWt::GetLongCommandPath(
 Wt::WTreeNode* item
)
{
  if (item == NULL) return "";

  // rebuild path:
  Wt::WString itemText = "";
  itemText = item->label()->text();

  while (item->parent() != NULL) {
    itemText = item->parentNode()->label()->text()+"/"+itemText;
    item = item->parentNode();
  }
  itemText = "/"+itemText;
  
  return itemText;
}

G4WTabWidget::G4WTabWidget(
):WTabWidget()
 ,tabSelected(false)
 ,lastCreated(-1)
{
}


  
void G4UIWt::TabCloseCallback(int a){
#if WT_VERSION >= 0x040500
  Wt::WWidget* temp = fTabWidget->widget(a);
  fTabWidget->closeTab (a);

  delete temp;

  if (fTabWidget->count() == 0) {
    if (fEmptyViewerTabLabel == NULL) {
      fEmptyViewerTabLabel = new Wt::WLabel("         If you want to have a Viewer, please use /vis/open commands. ");
    }

    fMyVSplitter->addWidget(fEmptyViewerTabLabel);
    fMyVSplitter->show();
    fEmptyViewerTabLabel->show();
    //    fTabWidget->setParent(0);
      fMainWindow->hide();
    delete fTabWidget;
    fTabWidget = NULL;
  }
#endif
}


void G4UIWt::HelpToolBoxActivated(){
  FillHelpTree();
}

// void G4UIWt::ToolBoxActivated(int a){

//   if (fToolBox->widget(a) == fHelpTBWidget) {
//     // Rebuild the help tree
//     FillHelpTree();
//   }
// }


Wt::WContainerWidget * G4UIWt::getLastTabContainerInsert() {
  return fLastTabWidgetInsert;
}


void G4WTabWidget::paintEvent(
Wt::WPaintDevice *
)
{
#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4WTabWidget::paintEvent OK\n");
#endif

  if (currentWidget()) {

    if ( isTabSelected()) {

      //      Wt::WCoreApplication::sendPostedEvents () ;

#ifdef G4DEBUG_INTERFACES_BASIC
      printf("G4WTabWidget::paintEvent OK\n");
#endif
      Wt::WString text = tabText (currentIndex()); 

      if (lastCreated == -1) {
        Wt::WString paramSelect = Wt::WString("/vis/viewer/select ")+text;
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if(UI != NULL)  {
          UI->ApplyCommand(paramSelect.toUTF8());
        }
      } else {
        lastCreated = -1;
      }
      setTabSelected(false);
    }
  }
}




#endif
