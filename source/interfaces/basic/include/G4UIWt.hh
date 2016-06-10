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
// $Id: G4UIWt.hh,v 1.23 2010-06-10 15:37:13 lgarnier Exp $
//
#ifndef G4UIWt_h
#define G4UIWt_h 

#if defined(G4UI_BUILD_WT_SESSION) || defined(G4UI_USE_WT)

#include <map>

#include "G4VBasicShell.hh"
#include "G4VInteractiveSession.hh"

#include <Wt/WObject>
#include <Wt/WWidget>
#include <Wt/WPushButton>
#include <Wt/WTree>
#include <Wt/WTreeNode>
#include <Wt/WMenu>
#include <Wt/WTabWidget>
#include <Wt/WStringListModel>



class G4UIsession;

// Class description :
//
//  G4UIWt : class to handle a Wt interactive session.
// G4UIWt is the Wt version of G4UIterminal.
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

class G4WTabWidget : public Wt::WTabWidget {
  public :
  G4WTabWidget();
  G4WTabWidget(Wt::WContainerWidget*&);
  inline void setTabSelected(bool a) { tabSelected = a; };
  inline void setLastTabCreated(int a) { lastCreated = a; };
  inline bool isTabSelected() { return tabSelected; };
  bool tabSelected;
  int lastCreated;
  int incTabPaint;

};


class G4UIWt : public Wt::WObject, public G4VBasicShell, public G4VInteractiveSession  {

public: // With description
  G4UIWt(int,char**);
  // (argv, argc) or (0, NULL) had to be given.
  G4UIsession* SessionStart();
  // To enter interactive Wt loop ; waiting/executing command,...
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

  bool AddTabWidget( Wt::WWidget*, Wt::WString,int,int);
  // To add a tab for vis openGL Qt driver
  
  Wt::WTabWidget* GetSceneTreeComponentsTBWidget();
  // Get the viewComponent

  bool IsSplitterReleased();

  inline bool IsIconMoveSelected() {
    return fMoveSelected;
  };
  inline bool IsIconRotateSelected() {
    return fRotateSelected;
  };
  inline bool IsIconPickSelected() {
    return fPickSelected;
  };
  inline bool IsIconZoomInSelected() {
    return fZoomInSelected;
  };
  inline bool IsIconZoomOutSelected() {
    return fZoomOutSelected;
  };

  /*  void SetIconMoveSelected();
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
   */

  inline  Wt::WContainerWidget * GetMainWindow() {
    return fMainWindow;
  };

public:
  ~G4UIWt();
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
  virtual void ExitHelp() const;

  void CreateHelpTree( Wt::WTreeNode*,G4UIcommandTree*);
  Wt::WTreeNode* FindTreeItem( Wt::WTreeNode *,const std::string&);

  Wt::WString GetCommandList(const G4UIcommand*);

  virtual G4bool GetHelpChoice(G4int&);// have to be implemeted because we heritate from G4VBasicShell
  bool eventFilter(Wt::WObject*,Wt::WEvent*);
  void ActivateCommand(G4String);
  //  QMap<int,Wt::WString> LookForHelpStringInChildTree(G4UIcommandTree *,const Wt::WString&);

  Wt::WContainerWidget* CreateVisParametersTBWidget();
  Wt::WWidget* CreateHelpTBWidget();
  Wt::WWidget* CreateCoutTBWidget();
  Wt::WWidget* CreateHistoryTBWidget();
  Wt::WWidget* CreateUITabWidget();
  Wt::WWidget* CreateSceneTreeComponentsTBWidget();
  Wt::WContainerWidget* CreateRightSplitterWidget();
  Wt::WContainerWidget* CreateLeftSplitterWidget();
  void OpenHelpTreeOnCommand(const Wt::WString &);
  Wt::WString GetShortCommandPath(const std::string & );
  Wt::WString GetLongCommandPath( Wt::WTreeNode*);
  G4bool IsGUICommand(const G4UIcommand*);
  bool CreateVisCommandGroupAndToolBox(G4UIcommand*, Wt::WWidget*, int, bool isDialog);
  bool CreateCommandWidget(G4UIcommand* command, Wt::WContainerWidget* parent, bool isDialog);

private:

  Wt::WContainerWidget * fMainWindow;
  Wt::WLabel *fCommandLabel;
  Wt::WLineEdit * fCommandArea;
  Wt::WTextArea *fCoutTBTextArea;
  Wt::WTextArea *fHelpArea;
  Wt::WTabWidget* fUITabWidget;
  Wt::WStringListModel fG4cout;
  Wt::WLineEdit * fCoutFilter;
  
  Wt::WSelectionBox *fHistoryTBTableList;
  Wt::WTree *fHelpTreeWidget;
  Wt::WPanel* fHelpTBWidget;
  Wt::WPanel* fHistoryTBWidget;
  Wt::WPanel* fCoutTBWidget;  
  Wt::WTabWidget* fSceneTreeComponentsTBWidget;
  Wt::WLineEdit* fHelpLine;
  G4WTabWidget* fViewerTabWidget;
  Wt::WString fCoutText;
  Wt::WLabel *fEmptyViewerTabLabel;
  Wt::WContainerWidget* fMainSplitterWidget;
  Wt::WContainerWidget* fRightSplitterWidget;
  Wt::WContainerWidget* fLeftSplitterWidget;
  Wt::WContainerWidget* fHelpVSplitter;
  
  Wt::WToolBar *fToolbarApp;
  Wt::WToolBar *fToolbarUser;
  Wt::WString fStringSeparator;
  G4String fLastErrMessage;
  Wt::WString fLastOpenPath;

  bool fMoveSelected;
  bool fRotateSelected;
  bool fPickSelected;
  bool fZoomInSelected;
  bool fZoomOutSelected;
  G4bool fExitSession;
  G4bool fExitPause;
  

  private :
  void ExitSession();
  void ClearButtonCallback();
  void CommandEnteredCallback();
  void CommandEditedCallback(const Wt::WString & text);
  void ButtonCallback(const char*);
  void HelpTreeClicCallback();
  void HelpTreeDoubleClicCallback();
  void ShowHelpCallback();
  void CommandHistoryCallback();
  void LookForHelpStringCallback();
  void CurrentChangedTabWidgetCallback(int);
  void CoutFilterCallback(const Wt::WString&);
  void TabCloseCallback(int);
  void ToolBoxActivated(int);
  void VisParameterCallback(Wt::WContainerWidget*);
  void ChangeColorCallback(Wt::WContainerWidget*);
  void ChangeCursorStyle(const Wt::WString&);
  void ChangeSurfaceStyle(const Wt::WString&);
  void OpenIconCallback(const Wt::WString&);
  void SaveIconCallback(const Wt::WString&);
  void ChangePerspectiveOrthoCallback(const Wt::WString&);
  
};

#endif

#endif

