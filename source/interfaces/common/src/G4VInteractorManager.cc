// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VInteractorManager.cc,v 1.2 1999-04-13 01:26:32 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

#include <stdlib.h>
#include <string.h>

#include "G4VInteractorManager.hh"

#define NewString(str)  \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : NULL)

/***************************************************************************/
G4VInteractorManager::G4VInteractorManager (
)
:argc(0)
,argv(NULL)
,mainInteractor(NULL)
,dispatchern(0)
,dispatchers(NULL)
,preActionn(0)
,preActions(NULL)
,postActionn(0)
,postActions(NULL)
,shelln(0)
,shells(NULL)
,secondaryLoopEnabled(TRUE)
,alreadyInSecondaryLoop(FALSE)
,exitSecondaryLoop(0)
,parentInteractor(NULL)
,createdInteractor(NULL)
,creationString(NULL)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
G4VInteractorManager::~G4VInteractorManager (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(argv!=NULL) {
    for(int argi=0;argi<argc;argi++) {
      if(argv[argi]!=NULL) free(argv[argi]);
    }
    free (argv);
  }
  argv                   = NULL;
  argc                   = 0;
  if(dispatchers!=NULL)  free(dispatchers);
  dispatchers            = NULL;
  dispatchern            = 0;
  if(preActions!=NULL)   free(preActions);
  preActions             = NULL;
  preActionn             = 0;
  if(postActions!=NULL)  free(postActions);
  postActions            = NULL;
  postActionn            = 0;
  if(shells!=NULL)       free(shells);
  shells                 = NULL;
  shelln                 = 0;

  secondaryLoopEnabled   = TRUE;
  alreadyInSecondaryLoop = FALSE;
  exitSecondaryLoop      = 0;
}
/***************************************************************************/
void G4VInteractorManager::SetArguments (
 int    a_argc
,char** a_argv
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  // Free previous values.
  if(argv!=NULL) {
    for(int argi=0;argi<argc;argi++) {
      if(argv[argi]!=NULL) free(argv[argi]);
    }
    free(argv);
  }
  argv = NULL;
  argc = 0;
  // Set new values.
  if(a_argc!=0) {
    argv = (char**)malloc(a_argc * sizeof(char*));
    if(argv!=NULL) {
      argc = a_argc;
      for(int argi=0;argi<a_argc;argi++) {
	argv[argi] = (char*)NewString (a_argv[argi]);
      }
    }
  }
}
/***************************************************************************/
char** G4VInteractorManager::GetArguments (
 int* a_argc
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_argc!=NULL) *a_argc = argc;
  return argv;
}
/***************************************************************************/
void G4VInteractorManager::SetMainInteractor (
 G4Interactor a_main
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  mainInteractor = a_main;
}
/***************************************************************************/
G4Interactor G4VInteractorManager::GetMainInteractor (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return mainInteractor;
}
/***************************************************************************/
void G4VInteractorManager::EnableSecondaryLoop (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  secondaryLoopEnabled = TRUE;
}
/***************************************************************************/
void G4VInteractorManager::DisableSecondaryLoop (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  secondaryLoopEnabled = FALSE;
}
/***************************************************************************/
void G4VInteractorManager::AddDispatcher (
 G4DispatchFunction a_dispatcher
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_dispatcher==NULL) return;
  for(int count=0;count<dispatchern;count++) {
    if(dispatchers[count]==a_dispatcher) return;  // Done.
  }
  if(dispatchern==0) dispatchers = (G4DispatchFunction*)malloc( sizeof(G4DispatchFunction));
  else               dispatchers = (G4DispatchFunction*)realloc(dispatchers, 
 				     (dispatchern+1) * sizeof(G4DispatchFunction));
  if(dispatchers==NULL) return;
  dispatchers[dispatchern] = a_dispatcher; 
  dispatchern++;
}
/***************************************************************************/
void G4VInteractorManager::RemoveDispatcher (
 G4DispatchFunction a_dispatcher
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  for(int count=0;count<dispatchern;count++) {
    if(dispatchers[count]==a_dispatcher) {
      dispatchers[count] = NULL;
      return;  // A dispatcher appears once in the List.
    }
  }
}
/***************************************************************************/
void G4VInteractorManager::DispatchEvent (
 void* a_event
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  for(int count=0;count<dispatchern;count++) {
    if( 
       (dispatchers[count]!=NULL)          &&
       (dispatchers[count](a_event)==TRUE) 
       )
      break;
  }
}
/***************************************************************************/
void G4VInteractorManager::AddSecondaryLoopPreAction (
 G4SecondaryLoopAction a_preAction
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_preAction==NULL) return;
  for(int count=0;count<preActionn;count++) {
    if(preActions[count]==a_preAction) return;  // Done.
  }
  if(preActionn==0) 
    preActions = (G4SecondaryLoopAction*)malloc( sizeof(G4SecondaryLoopAction));
  else
    preActions = (G4SecondaryLoopAction*)realloc(preActions, 
				(preActionn+1) * sizeof(G4SecondaryLoopAction));
  if(preActions==NULL) return;
  preActions[preActionn] = a_preAction; 
  preActionn++;
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoopPreActions (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  for(int count=0;count<preActionn;count++) {
    if(preActions[count]!=NULL) preActions[count]();
  }
}
/***************************************************************************/
void G4VInteractorManager::AddSecondaryLoopPostAction (
 G4SecondaryLoopAction a_postAction
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_postAction==NULL) return;
  for(int count=0;count<postActionn;count++) {
    if(postActions[count]==a_postAction) return;  // Done.
  }
  if(postActionn==0) 
    postActions = (G4SecondaryLoopAction*)malloc( sizeof(G4SecondaryLoopAction));
  else
    postActions = (G4SecondaryLoopAction*)realloc(postActions, 
   			         (postActionn+1) * sizeof(G4SecondaryLoopAction));
  if(postActions==NULL) return;
  postActions[postActionn] = a_postAction; 
  postActionn++;
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoopPostActions (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  for(int count=0;count<postActionn;count++) {
    if(postActions[count]!=NULL) postActions[count]();
  }
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoop (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(Inited()==FALSE) return;

  if(secondaryLoopEnabled==FALSE) return;
  
  if (alreadyInSecondaryLoop==FALSE) {
    alreadyInSecondaryLoop   = TRUE;
    exitSecondaryLoop        = 0;
    SecondaryLoopPreActions  ();
    //      for(int count=0;count<shelln;count++) XWidgetUniconify(shells[count]);
    void*                    event;
    while(1) {
      event = GetEvent();
      if(event==NULL) break;
      DispatchEvent  (event);
      if(exitSecondaryLoop!=0) break;
    }
    SecondaryLoopPostActions ();
    }
}
/***************************************************************************/
void G4VInteractorManager::RequireExitSecondaryLoop (
 int a_code
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(secondaryLoopEnabled==FALSE) return;
  if(a_code==0)            a_code = 1;
  exitSecondaryLoop        = a_code;
  alreadyInSecondaryLoop   = FALSE;
  // for(int count=0;count<shelln;count++) XWidgetIconify(shells[count]);
  // if(shelln!=0)            XSync(XtDisplay(topWidget),False);
}
/***************************************************************************/
int G4VInteractorManager::GetExitSecondaryLoopCode (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return exitSecondaryLoop;
}
/***************************************************************************/
void G4VInteractorManager::AddShell (
 G4Interactor a_shell
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_shell==NULL) return;
  for(int count=0;count<shelln;count++) {
    if(shells[count]==a_shell) return;  // Done.
  }
  if(shelln==0)    shells = (G4Interactor*)malloc( sizeof(G4Interactor));
  else             shells = (G4Interactor*)realloc(shells, (shelln+1) * sizeof(G4Interactor));
  if(shells==NULL) return;
  shells[shelln] = a_shell; 
  shelln++;
}
/***************************************************************************/
void G4VInteractorManager::RemoveShell (
 G4Interactor a_shell
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  for(int count=0;count<shelln;count++) {
    if(shells[count]==a_shell) {
      shells[count] = NULL;
      return;  // A shell appears once in the List.
    }
  }
}
/***************************************************************************/
void G4VInteractorManager::SetParentInteractor (
 G4Interactor a_interactor
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  parentInteractor = a_interactor;
}
/***************************************************************************/
G4Interactor G4VInteractorManager::GetParentInteractor (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return parentInteractor;
}
/***************************************************************************/
void G4VInteractorManager::SetCreatedInteractor (
 G4Interactor a_interactor
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  createdInteractor = a_interactor;
}
/***************************************************************************/
G4Interactor G4VInteractorManager::GetCreatedInteractor (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return createdInteractor;
}
/***************************************************************************/
void G4VInteractorManager::SetCreationString (
 char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  creationString = a_string;
}
/***************************************************************************/
char* G4VInteractorManager::GetCreationString (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return creationString;
}

