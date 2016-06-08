// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include <string.h>
#include <stdlib.h>

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4VInteractiveSession.hh"

#include "G4InteractorMessenger.hh"

#define STRDUP(str)  ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : (char*)NULL)
#define STRDEL(str) {if((str)!=NULL) {free(str);str=NULL;}}

static G4bool GetValues (G4String,int,G4String*); 

G4InteractorMessenger::G4InteractorMessenger (
 G4VInteractiveSession* a_session
)
{
  session = a_session;

  G4UIparameter* parameter;

  interactorDirectory = new G4UIdirectory("/interactor/");
  interactorDirectory->SetGuidance("UI interactors commands.");

  // /interactor/addMenu :
  addMenu = new G4UIcommand("/interactor/addMenu",this);
  addMenu->SetGuidance("Add a menu to menu bar.");
  parameter = new G4UIparameter("Name",'s',false);
  parameter->SetDefaultValue("dummy");
  addMenu->SetParameter (parameter);
  parameter = new G4UIparameter("Label",'s',false);
  parameter->SetDefaultValue("dummy");
  addMenu->SetParameter (parameter);

  // /interactor/addButton :
  addButton = new G4UIcommand("/interactor/addButton",this);
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

}

G4InteractorMessenger::~G4InteractorMessenger()
{
  delete addButton;
  delete addMenu;
  delete interactorDirectory;
}

void G4InteractorMessenger::SetNewValue (
 G4UIcommand* command
,G4String newValue
) 
{
  int paramn = command->GetParameterEntries();
  G4String* params = new G4String [paramn];
  if(GetValues(newValue,paramn,params)==true) {
    if(command==addMenu) {
      session->AddMenu((const char*)params[0],(const char*)params[1]);
    } else if(command==addButton) {
      session->AddButton((const char*)params[0],(const char*)params[1],(const char*)params[2]);
    }
  }
  delete [] params;
}
G4bool GetValues (
 G4String newValue
,int paramn
,G4String* params
) 
{
  char* value = STRDUP(newValue.data());
  if(value==NULL) return false;
  char* tok = strtok(value," ");
  for( int i=0; i<paramn;i++ ) {
    if(tok==NULL) {
      STRDEL(value);
      return false;
    }
    G4String token = tok;
    if( token(0)=='"' ) {
      while( token(token.length()-1) != '"' ) {
	tok = strtok(NULL," ");
	if( (tok==NULL) || (*tok=='\0')) {
	  STRDEL(value);
	  return false;
	}
	token += " ";
	token += tok;
      }
      token = (G4String)token.strip(G4String::both,'"');
    }
    if( token.isNull() ) {
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


