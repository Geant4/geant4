// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VInteractorManager.cc,v 1.4 1999-05-11 12:33:54 barrand Exp $
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
  argv = NULL;
  argc = 0;
  dispatchers.clear();
  preActions.clear();
  postActions.clear();
  shells.clear();
  secondaryLoopEnabled = TRUE;
  alreadyInSecondaryLoop = FALSE;
  exitSecondaryLoop = 0;
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
  if(dispatchers.contains(a_dispatcher)==TRUE) return;
  dispatchers.insert(a_dispatcher);
}
/***************************************************************************/
void G4VInteractorManager::RemoveDispatcher (
 G4DispatchFunction a_dispatcher
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  dispatchers.removeAll(a_dispatcher);
}
/***************************************************************************/
void G4VInteractorManager::DispatchEvent (
 void* a_event
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  int dispatchern = dispatchers.entries();
  G4DispatchFunction func;
  G4bool status;
  for(int count=0;count<dispatchern;count++) {
    func = dispatchers[count];
    if(func!=NULL) {
      if(func(a_event)==true) return;
    }
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
  if(preActions.contains(a_preAction)==TRUE) return;
  preActions.insert(a_preAction);
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoopPreActions (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  int preActionn = preActions.entries();
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
  if(postActions.contains(a_postAction)==TRUE) return;
  postActions.insert(a_postAction);
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoopPostActions (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  int postActionn = postActions.entries();
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
    //for(int count=0;count<shelln;count++) XWidgetUniconify(shells[count]);
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
  if(shells.contains(a_shell)==TRUE) return;
  shells.insert(a_shell);
}
/***************************************************************************/
void G4VInteractorManager::RemoveShell (
 G4Interactor a_shell
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  shells.removeAll(a_shell);
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

