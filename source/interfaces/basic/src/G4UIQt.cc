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
// $Id: G4UIQt.cc,v 1.4 2007-11-09 15:03:21 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// L. Garnier

//#define DEBUG

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
#include <qmainwindow.h>
#include <qlineedit.h>
#include <qwidget.h>
#include <qmenu.h>
#include <qmenubar.h>
#include <qboxlayout.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qsplitter.h>
#include <qscrollbar.h>
#include <qdialog.h>
#include <qevent.h>
#include <qlistwidget.h>
#include <qtextedit.h>
#include <qtreewidget.h>
#include <qsignalmapper.h>


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
  :fHelpDialog(NULL)
{
  G4Qt* interactorManager = G4Qt::getInstance ();
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  fMainWindow = new QMainWindow();
  fMainWindow->setWindowTitle( "G4UI Session" ); 
  fMainWindow->resize(800,600); 
  fMainWindow->move(QPoint(50,100));

  QSplitter *splitter = new QSplitter(Qt::Vertical);
  fTextArea = new QTextEdit();
  QPushButton *clearButton = new QPushButton("clear");
  connect(clearButton, SIGNAL(clicked()), SLOT(ClearButtonCallback()));

  fCommandHistoryArea = new QListWidget();
  fCommandHistoryArea->setSelectionMode(QAbstractItemView::SingleSelection);
  connect(fCommandHistoryArea, SIGNAL(itemSelectionChanged()), SLOT(CommandHistoryCallback()));
  fCommandHistoryArea->installEventFilter(this);
  fCommandLabel = new QLabel();

  fCommandArea = new QLineEdit();
  fCommandArea->installEventFilter(this);
  fCommandArea->activateWindow();
  connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));
  fCommandArea->setFocusPolicy ( Qt::StrongFocus );
  fCommandArea->setFocus(Qt::TabFocusReason);
  fTextArea->setReadOnly(true);


  // Set layouts

  QWidget* topWidget = new QWidget();
  QVBoxLayout *layoutTop = new QVBoxLayout;

  QWidget* bottomWidget = new QWidget();
  QVBoxLayout *layoutBottom = new QVBoxLayout;


  layoutTop->addWidget(fTextArea);
  layoutTop->addWidget(clearButton);
  topWidget->setLayout(layoutTop);

  layoutBottom->addWidget(fCommandHistoryArea);
  layoutBottom->addWidget(fCommandLabel);
  layoutBottom->addWidget(fCommandArea);
  bottomWidget->setLayout(layoutBottom);


  splitter->addWidget(topWidget);
  splitter->addWidget(bottomWidget);
  fMainWindow->setCentralWidget(splitter);

  // Add a quit subMenu
  QMenu *fileMenu = fMainWindow->menuBar()->addMenu("File");
  fileMenu->addAction("Quitter", fMainWindow, SLOT(close()));

  // Add a Help menu
  QMenu *helpMenu = fMainWindow->menuBar()->addMenu("Help");
  helpMenu->addAction("Show Help", this, SLOT(ShowHelpCallback()));

  // Set the splitter size. The fTextArea sould be 2/3 on the fMainWindow
  QList<int> vals = splitter->sizes();
  if(vals.size()==2) {
    vals[0] = (splitter->orientation()==Qt::Vertical ? splitter->height() : splitter->width())*3/4;
    vals[1] = (splitter->orientation()==Qt::Vertical ? splitter->height() : splitter->width())*1/4;
    splitter->setSizes(vals);
  }


  if(UI!=NULL) UI->SetCoutDestination(this);  // TO KEEP
}



G4UIQt::~G4UIQt(
) 
{ 
  G4UImanager* UI = G4UImanager::GetUIpointer();  // TO KEEP
  if(UI!=NULL) {  // TO KEEP
    UI->SetSession(NULL);  // TO KEEP
    UI->SetCoutDestination(NULL);  // TO KEEP
  }
  
  if (fMainWindow!=NULL)
    delete fMainWindow;
}



/**   Start the Qt main loop
*/
G4UIsession* G4UIQt::SessionStart (
)
{

  G4Qt* interactorManager = G4Qt::getInstance ();
  fMainWindow->show();
  Prompt("session");
  exitSession = false;


  printf("disable secondary loop\n");
  interactorManager->DisableSecondaryLoop (); // TO KEEP
  ((QApplication*)interactorManager->GetMainInteractor())->exec(); 
  // on ne passe pas le dessous ? FIXME ????
  // je ne pense pas 13/06

  //   void* event; // TO KEEP
  //   while((event = interactorManager->GetEvent())!=NULL) {  // TO KEEP
  //     interactorManager->DispatchEvent(event); // TO KEEP
  //     if(exitSession==true) break; // TO KEEP
  //   } // TO KEEP

  interactorManager->EnableSecondaryLoop ();
  printf("enable secondary loop\n");
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

  printf("G4UIQt::PauseSessionStart\n");
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

  printf("G4UIQt::SecondaryLoop\n");
  G4Qt* interactorManager = G4Qt::getInstance (); // TO KEEP ?
  Prompt(aPrompt); // TO KEEP
  exitPause = false; // TO KEEP
  void* event; // TO KEEP
  while((event = interactorManager->GetEvent())!=NULL) {  // TO KEEP
    interactorManager->DispatchEvent(event); // TO KEEP
    if(exitPause==true) break; // TO KEEP
  } // TO KEEP
  Prompt("session"); // TO KEEP
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
  G4Qt* interactorManager = G4Qt::getInstance ();
  if (!interactorManager) return 0;
  
  //  printf(" **************** G4 Cout : %s\n",(char*)aString.data());
  fTextArea->append(QString((char*)aString.data()).trimmed());
  fTextArea->verticalScrollBar()->setSliderPosition(fTextArea->verticalScrollBar()->maximum());
  interactorManager->FlushAndWaitExecution();
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
  G4Qt* interactorManager = G4Qt::getInstance ();
  if (!interactorManager) return 0;

  QColor previousColor = fTextArea->textColor();
  fTextArea->setTextColor(Qt::red);
  fTextArea->append(QString((char*)aString.data()).trimmed());
  fTextArea->setTextColor(previousColor);
  fTextArea->verticalScrollBar()->setSliderPosition(fTextArea->verticalScrollBar()->maximum());
  interactorManager->FlushAndWaitExecution();
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
  signalMapper->setMapping(action, QString(aCommand));
  connect(signalMapper, SIGNAL(mapped(const QString &)),this, SLOT(ButtonCallback(const QString&)));
}




/**
   Open the help dialog in a separate window.<br>
   This will be display as a tree widget.<br>
   Implementation of <b>void G4VBasicShell::TerminalHelp(G4String newCommand)</b>

   @param newCommand : open the tree widget item on this command if is set
*/
void G4UIQt::TerminalHelp(
 G4String newCommand
)
{
  // Create the help dialog
  if (!fHelpDialog) {
    fHelpDialog = new QDialog(fMainWindow,Qt::WindowTitleHint | Qt::WindowSystemMenuHint | Qt::WindowMinMaxButtonsHint);
    QSplitter *splitter = new QSplitter(Qt::Horizontal);
    fHelpArea = new QTextEdit();
    QPushButton *exitButton = new QPushButton("Exit");
    connect(exitButton, SIGNAL(clicked()), fHelpDialog,SLOT(close()));
    fHelpArea->setReadOnly(true);

    // the help tree
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if(UI==NULL) return;
    G4UIcommandTree * treeTop = UI->GetTree();

    // build widget
    fHelpTreeWidget = new QTreeWidget();
    fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    fHelpTreeWidget->setColumnCount(2);
    fHelpTreeWidget->setColumnHidden(1,true);
    QStringList labels;
    labels << QString("Command") << QString("Description");
    fHelpTreeWidget->setHeaderLabels(labels);

    QList<QTreeWidgetItem *> items;
    G4int treeSize = treeTop->GetTreeEntry();
    QTreeWidgetItem * newItem;
    for (int a=0;a<treeSize;a++) {
      // Creating new item
      QStringList stringList;
      stringList << QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).trimmed()  ;
      stringList << QString((char*)(treeTop->GetTree(a+1)->GetTitle()).data()).trimmed()  ;

      newItem = new QTreeWidgetItem(stringList);

      // look for childs
      CreateChildTree(newItem,treeTop->GetTree(a+1));
      items.append(newItem);
    }
    fHelpTreeWidget->insertTopLevelItems(0, items);

    connect(fHelpTreeWidget, SIGNAL(itemSelectionChanged ()),this, SLOT(HelpTreeClicCallback()));  
    connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked (QTreeWidgetItem*,int)),this, SLOT(HelpTreeDoubleClicCallback(QTreeWidgetItem*,int)));  

    // Set layouts

    QVBoxLayout *vLayout = new QVBoxLayout;

    splitter->addWidget(fHelpTreeWidget);
    splitter->addWidget(fHelpArea);

    vLayout->addWidget(splitter);
    vLayout->addWidget(exitButton);
    fHelpDialog->setLayout(vLayout);

  }

  // Look for the choosen command "newCommand"
  size_t i = newCommand.index(" ");
  G4String targetCom="";
  if( i != std::string::npos )
    {
      G4String newValue = newCommand(i+1,newCommand.length()-(i+1));
      newValue.strip(G4String::both);
      targetCom = ModifyToFullPathCommand( newValue );
    }
  if (targetCom != "") {
    QTreeWidgetItem* findItem = NULL;
    for (int a=0;a<fHelpTreeWidget->topLevelItemCount();a++) {
      if (!findItem) {
        findItem = FindTreeItem(fHelpTreeWidget->topLevelItem(a),QString((char*)targetCom.data()));
      }
    }
    
    if (findItem) {      
      
      //collapsed open item
      QList<QTreeWidgetItem *> selected;
      selected = fHelpTreeWidget->selectedItems();
      if ( selected.count() != 0 ) {
        QTreeWidgetItem * tmp =selected.at( 0 );
        while ( tmp) {
          tmp->setExpanded(false);
          tmp = tmp->parent();
        }
      }
      
      // clear old selection
      fHelpTreeWidget->clearSelection();

      // set new selection
      findItem->setSelected(true);
      
      // expand parent item
      while ( findItem) {
        findItem->setExpanded(true);
        findItem = findItem->parent();
      }

      // Call the update of the right textArea
      HelpTreeClicCallback();
    }
  }
  fHelpDialog->setWindowTitle("Help on commands"); 
  fHelpDialog->resize(800,600); 
  fHelpDialog->move(QPoint(400,150));
  fHelpDialog->show();
  fHelpDialog->raise();
  fHelpDialog->activateWindow();
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

  // Get the Sub directories
  for (int a=0;a<aCommandTree->GetTreeEntry();a++) {

    QStringList stringList;
    stringList << QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).trimmed()  ;
    stringList << QString((char*)(aCommandTree->GetTree(a+1)->GetTitle()).data()).trimmed()  ;
    newItem = new QTreeWidgetItem(stringList);

    CreateChildTree(newItem,aCommandTree->GetTree(a+1));
    aParent->addChild(newItem);
  }



  // Get the Commands

  for (int a=0;a<aCommandTree->GetCommandEntry();a++) {
    
    QStringList stringList;
    stringList << QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed()  ;
    stringList << QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed()  ;
    newItem = new QTreeWidgetItem(stringList);
      
    aParent->addChild(newItem);
    newItem->setExpanded(false);
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

  if (aParent->text(0) == aCommand)
    return aParent;
  
  QTreeWidgetItem * tmp = NULL;
  for (int a=0;a<aParent->childCount();a++) {
    if (!tmp)
      tmp = FindTreeItem(aParent->child(a),aCommand);
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
      txt += " Parameter type  : " + QString(param->GetParameterType())+ "\n";
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
 G4int& aInt
)
{
  printf("G4UIQt::GetHelpChoice SHOULD NEVER GO HERE");
  return true;
}


/**   Implement G4VBasicShell vurtual function
*/
void G4UIQt::ExitHelp(
)
{
  printf("G4UIQt::ExitHelp SHOULD NEVER GO HERE");
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
  if (aObj == NULL) return false;
  if (aEvent == NULL) return false;

  if (aObj == fCommandHistoryArea) {
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
        int selection = fCommandHistoryArea->currentRow();
        if (fCommandHistoryArea->count()) {
          if (selection == -1) {
            selection = fCommandHistoryArea->count()-1;
          } else {
            if (e->key() == (Qt::Key_Down)) {
              if (selection <(fCommandHistoryArea->count()-1))
                selection++;
            } else if (e->key() == (Qt::Key_PageDown)) {
              selection = fCommandHistoryArea->count()-1;
            } else if (e->key() == (Qt::Key_Up)) {
              if (selection >0)
                selection --;
            } else if (e->key() == (Qt::Key_PageUp)) {
              selection = 0;
            }
          }
          fCommandHistoryArea->clearSelection();
          fCommandHistoryArea->item(selection)->setSelected(true);
          fCommandHistoryArea->setCurrentItem(fCommandHistoryArea->item(selection));
        }
      } else if (e->key() == (Qt::Key_Tab)) {
        G4String ss = Complete(fCommandArea->text().toStdString().c_str());
        fCommandArea->setText((char*)(ss.data()));

        // do not pass by parent, it will disable widget tab focus !
        return true;
      }
    }
  }
  // pass the event on to the parent class
  return QObject::eventFilter(aObj, aEvent);
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
  fTextArea->clear();
}


/**   Callback call when "click on a menu entry.<br>
   Send the associated command to geant4
*/
void G4UIQt::CommandEnteredCallback (
)
{
  G4String command (fCommandArea->text().toStdString().c_str());
  if (fCommandArea->text().trimmed() != "") {
    fCommandHistoryArea->addItem(fCommandArea->text());
    fCommandHistoryArea->clearSelection();
    fCommandHistoryArea->setCurrentItem(NULL);
    fCommandArea->setText("");

    G4Qt* interactorManager = G4Qt::getInstance ();
    if (interactorManager) { 
      interactorManager->FlushAndWaitExecution();
    }
    if (command(0,4) != "help") {
      ApplyShellCommand (command,exitSession,exitPause);
    } else {
      TerminalHelp(command);
    }
    printf("after \n");
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
  if(exitSession==true) 
    SessionTerminate();
}



/**   This callback is activated when user selected a item in the help tree
*/
void G4UIQt::HelpTreeClicCallback (
)
{
  //  printf("G4UIQt::HelpTreeClicCallback");
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
  G4UIcommand* command = treeTop->FindPath(item->text (1).toStdString().c_str());
  if (command) {
    fHelpArea->setText(GetCommandList(command));
  } else {
    // this is not a command, this is a sub directory
    // We display the Title
    fHelpArea->setText(item->text (1).toStdString().c_str());
  }
}

/**   This callback is activated when user double clic on a item in the help tree
*/
void G4UIQt::HelpTreeDoubleClicCallback (
QTreeWidgetItem* item
 ,int nb
)
{
  printf("G4UIQt::HelpTreeDoubleClicCallback");
  HelpTreeClicCallback();
  fCommandArea->setText(item->text (1));
}


/**   Callback called when user select an old command in the command history<br>
   Give it to the command area.
*/
void G4UIQt::CommandHistoryCallback(
)
{
  QListWidgetItem* item =  NULL;
  if (!fCommandHistoryArea)
    return ;

  
  QList<QListWidgetItem *> list =fCommandHistoryArea->selectedItems();
  if (list.isEmpty())
    return;
  item = list.first();
  if (!item)
    return;
  fCommandArea->setText(item->text());

}

#endif
