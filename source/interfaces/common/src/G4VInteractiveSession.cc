// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4InteractorMessenger.hh"

#include "G4VInteractiveSession.hh"

/***************************************************************************/
G4VInteractiveSession::G4VInteractiveSession (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  messenger = new G4InteractorMessenger(this);
}
/***************************************************************************/
G4VInteractiveSession::~G4VInteractiveSession() 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  delete messenger;
}
/***************************************************************************/
void G4VInteractiveSession::AddMenu (
 const char* a_name
,const char* a_label
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4VInteractiveSession::AddButton (
 const char* a_menu
,const char* a_label
,const char* a_command
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4VInteractiveSession::AddInteractor (
 G4String a_name
,G4Interactor a_interactor
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  interactors[a_name] = a_interactor;
}
/***************************************************************************/
G4Interactor G4VInteractiveSession::GetInteractor (
 G4String a_name
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4interactor_map::iterator it;
  if((it=interactors.find(a_name))==interactors.end()) return NULL;
  return (*it).second;  
}


