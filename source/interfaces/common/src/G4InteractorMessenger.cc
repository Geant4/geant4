// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include <rw/ctoken.h>

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4VInteractiveSession.hh"

#include "G4InteractorMessenger.hh"

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
  RWCTokenizer newValueToken( newValue );
  G4String aToken;
  for( int i=0; i<paramn;i++ ) {
    aToken = (G4String)newValueToken();
    if( aToken(0)=='"' ) {
      while( aToken(aToken.length()-1) != '"' ) {
	G4String additionalToken = (G4String)newValueToken();
	if( additionalToken.isNull() ) { 
	  return false;
	}
	aToken += " ";
	aToken += additionalToken;
      }
      aToken = (G4String)aToken.strip(G4String::both,'"');
    }
    if( aToken.isNull() ) {
      return false;
    } else { 
      params[i] = aToken;
    }
  }
  return true;
}


