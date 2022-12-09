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

#include <string.h>
#include <stdlib.h>

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4VInteractiveSession.hh"

#include "G4InteractorMessenger.hh"

#define STRDUP(str)  ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : (char*)NULL)
#define STRDEL(str) {if((str)!=NULL) {free(str);str=NULL;}}

static G4bool GetValues (G4String, G4int, G4String*); 

G4InteractorMessenger::G4InteractorMessenger (
 G4VInteractiveSession* a_session
)
{
  session = a_session;

  G4UIparameter* parameter;

  // gui commands should *not* be broadcast to workers
  G4bool propagateToWorkers;
  interactorDirectory = new G4UIdirectory("/gui/",propagateToWorkers=false);
  interactorDirectory->SetGuidance("UI interactors commands.");

  // /gui/addMenu :
  addMenu = new G4UIcommand("/gui/addMenu",this);
  addMenu->SetGuidance("Add a menu to menu bar.");
  parameter = new G4UIparameter("Name",'s',false);
  parameter->SetDefaultValue("dummy");
  addMenu->SetParameter (parameter);
  parameter = new G4UIparameter("Label",'s',false);
  parameter->SetDefaultValue("dummy");
  addMenu->SetParameter (parameter);

  // /gui/addButton :
  addButton = new G4UIcommand("/gui/addButton",this);
  addButton->SetGuidance("Add a button to menu.");
  parameter = new G4UIparameter("Menu",'s',false);
  parameter->SetDefaultValue("dummy");
  addButton->SetParameter (parameter);
  parameter = new G4UIparameter("Label",'s',false);
  parameter->SetDefaultValue("dummy");
  addButton->SetParameter (parameter);
  parameter = new G4UIparameter("Command",'s',false);
  parameter->SetDefaultValue("");
  addButton->SetParameter (parameter);

  // /gui/defaultIcons :
  defaultIcons = new G4UIcommand("/gui/defaultIcons",this);
  defaultIcons->SetGuidance("Set the Geant4 defaults icons in Qt driver.");
  defaultIcons->SetGuidance("By default, Geant4 icons are enable.");

  parameter = new G4UIparameter("bool",'b',true);
  parameter->SetDefaultValue("true");
  defaultIcons->SetParameter (parameter);
  
  // /gui/addIcon :
  addIcon = new G4UIcommand("/gui/addIcon",this);
  addIcon->SetGuidance
  ("Add a non-checkable icon to the Icon toolbar.");
  addIcon->SetGuidance
  ("If the Icon parameter is set to \"user_icon\", you should provide the icon file in xpm format, otherwise you have to choose one of the candidate icons");
  addIcon->SetGuidance
  ("A command given without parameters will display a window that will allow one to choose the parameters (if needed) for this command.");
  addIcon->SetGuidance
  ("E.g: /gui/addIcon \"Change background color\" user_icon /vis/viewer/set/background ../Images/background.xpm");
  addIcon->SetGuidance
  ("Special cases for the Icon parameter:");
  addIcon->SetGuidance
  (" - open: Open an open-file-selector that can run the Command with File as argument.");
  addIcon->SetGuidance
  (" - save: Open a save-file-selector that can run the Command with File as argument.");
  addIcon->SetGuidance
  (" - move/rotate/pick/zoom_in/zoom_out: Theses icons are radio-button icons that can change cursor action.");
  addIcon->SetGuidance
  (" - wireframe/solid/hidden_line_removal/hidden_line_and_surface_removal: These icons are radio-button icons that can change drawing style.");
  addIcon->SetGuidance
  (" - perspective/ortho: These icons are radio-button icons that can change projection style.");

  parameter = new G4UIparameter("Label",'s',false);
  parameter->SetDefaultValue("");
  addIcon->SetParameter (parameter);

  parameter = new G4UIparameter("Icon",'s',false);
  parameter->SetDefaultValue("");
  parameter->SetParameterCandidates
  ("open save move rotate pick zoom_in zoom_out wireframe solid hidden_line_removal hidden_line_and_surface_removal perspective ortho exit user_icon");
  addIcon->SetParameter (parameter);

  parameter = new G4UIparameter("Command",'s',true);
  parameter->SetDefaultValue("no_command");
  addIcon->SetParameter (parameter);

  parameter = new G4UIparameter("File",'s',true);
  parameter->SetDefaultValue("no_file");
  addIcon->SetParameter (parameter);

  // /gui/system :
  sys = new G4UIcommand("/gui/system",this);
  sys->SetGuidance("Send a command to the system.");
  parameter = new G4UIparameter("Command",'s',false);
  parameter->SetDefaultValue("");
  sys->SetParameter (parameter);

  // /gui/outputStyle :
  outputStyle = new G4UIcommand("/gui/outputStyle",this);
  outputStyle->SetGuidance("Set output style.");
  outputStyle->SetGuidance("Highlights commands if requested and if /control/verbose > 0.");
  parameter = new G4UIparameter("destination",'s',true);  // Omitable
  parameter->SetParameterCandidates("cout cerr warnings errors all");
  parameter->SetDefaultValue("all");
  outputStyle->SetParameter (parameter);
  parameter = new G4UIparameter("type",'s',true);  // Omitable
  parameter->SetParameterCandidates("fixed proportional");
  parameter->SetDefaultValue("fixed");
  outputStyle->SetParameter (parameter);
  parameter = new G4UIparameter("highlight",'s',true);  // Omitable
  parameter->SetParameterCandidates("highlight no-highlight");
  parameter->SetDefaultValue("highlight");
  outputStyle->SetParameter (parameter);

  // /gui/nativeMenuBar :
  nativeMenu = new G4UIcommand("/gui/nativeMenuBar",this);
  nativeMenu->SetGuidance("Allow native menu bar in Geant4 Qt driver.");
  nativeMenu->SetGuidance("By default, enable.");

  parameter = new G4UIparameter("bool",'b',true);
  parameter->SetDefaultValue("true");
  nativeMenu->SetParameter (parameter);
  // /gui/clearMenu
  clearMenu = new G4UIcommand("/gui/clearMenu",this);
  clearMenu->SetGuidance("Clear menu bar, remove all user defined menu entries.");
}

G4InteractorMessenger::~G4InteractorMessenger()
{
  delete clearMenu;
  delete nativeMenu;
  delete outputStyle;
  delete sys;
  delete defaultIcons;
  delete addIcon;
  delete addButton;
  delete addMenu;
  delete interactorDirectory;
}

void G4InteractorMessenger::SetNewValue (
 G4UIcommand* command
,G4String newValue
) 
{
  G4int paramn = (G4int)command->GetParameterEntries();
  G4String* params = new G4String [paramn];
  if(GetValues(newValue,paramn,params)==true) {
    if(command==addMenu) {
      session->AddMenu((const char*)params[0],(const char*)params[1]);
    } else if(command==addButton) {
      session->AddButton((const char*)params[0],(const char*)params[1],(const char*)params[2]);
    } else if(command==addIcon) {
      session->AddIcon((const char*)params[0],(const char*)params[1],(const char*)params[2],(const char*)params[3]);
    } else if(command==defaultIcons) {
      session->DefaultIcons(command->ConvertToBool(newValue));
    } else if(command==sys) {
      G4int rc = system((const char*)params[0]);
      if ( rc < 0 ){ } 
    } else if(command==outputStyle) {
      session->OutputStyle((const char*)params[0],(const char*)params[1],(const char*)params[2]);
    } else if(command==nativeMenu) {
      session->NativeMenu(command->ConvertToBool(newValue));
    } else if(command==clearMenu) {
      session->ClearMenu();
    }
  }
  delete [] params;
}
G4bool GetValues (
 G4String newValue
,G4int paramn
,G4String* params
) 
{
  char* value = STRDUP(newValue.data());
  if(value==NULL) return false;
  char* tok = strtok(value," ");
  for( G4int i=0; i<paramn; ++i ) {
    if(tok==NULL) {
      STRDEL(value);
      return false;
    }
    G4String token = tok;
    if( token[0]=='"' ) {
      while( token.back() != '"' ) {
	tok = strtok(NULL," ");
	if( (tok==NULL) || (*tok=='\0')) {
	  STRDEL(value);
	  return false;
	}
	token += " ";
	token += tok;
      }
      G4StrUtil::strip(token, '"');
    }
    if( token.empty() ) {
      STRDEL(value);
      return false;
    } else { 
      params[i] = token;
    }
    tok = strtok(NULL," ");
  }
  STRDEL(value);
  return true;
}
