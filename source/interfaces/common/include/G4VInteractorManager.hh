// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VInteractorManager.hh,v 1.2 1999-04-13 01:26:30 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

#ifndef G4VINTERACTORMANAGER_HH
#define G4VINTERACTORMANAGER_HH

#include "globals.hh"

typedef void*  G4Interactor;
typedef G4bool (*G4DispatchFunction)(void*);
typedef void   (*G4SecondaryLoopAction)();

class G4VInteractorManager {
public:
                 G4VInteractorManager            ();                     
  virtual       ~G4VInteractorManager            ();                     
  void           SetArguments                    (int,char**);
  char**         GetArguments                    (int*);
  void           SetMainInteractor               (G4Interactor);
  G4Interactor   GetMainInteractor               ();
  void           AddDispatcher                   (G4DispatchFunction);
  void           RemoveDispatcher                (G4DispatchFunction);
  void           AddSecondaryLoopPreAction       (G4SecondaryLoopAction);
  void           AddSecondaryLoopPostAction      (G4SecondaryLoopAction);
  void           AddShell                        (G4Interactor);
  void           RemoveShell                     (G4Interactor);
  void           EnableSecondaryLoop             ();
  void           DisableSecondaryLoop            ();
  void           SecondaryLoopPreActions         ();
  void           SecondaryLoopPostActions        ();
  void           RequireExitSecondaryLoop        (int); 
  void           DispatchEvent                   (void*);
  void           SecondaryLoop                   (); 
  int            GetExitSecondaryLoopCode        ();
  void           PutStringInResourceDatabase     (char*);
  virtual G4bool Inited                          () = 0;
  virtual void*  GetEvent                        () = 0;
  virtual void   FlushAndWaitExecution           () = 0;
  void           SetParentInteractor             (G4Interactor);
  G4Interactor   GetParentInteractor             ();
  void           SetCreatedInteractor            (G4Interactor);
  G4Interactor   GetCreatedInteractor            ();
  void           SetCreationString               (char*);
  char*          GetCreationString               ();
private:
  int                    argc;
  char**                 argv;
  G4Interactor           mainInteractor;
  int                    dispatchern;
  G4DispatchFunction*    dispatchers;
  int                    preActionn;
  G4SecondaryLoopAction* preActions;
  int                    postActionn;
  G4SecondaryLoopAction* postActions;
  int                    shelln;
  G4Interactor*          shells;
  G4bool                 secondaryLoopEnabled;
  G4bool                 alreadyInSecondaryLoop;
  int                    exitSecondaryLoop;
  G4Interactor           parentInteractor;
  G4Interactor           createdInteractor;
  char*                  creationString;
};

#define OGL_EXIT_CODE 1
#define OIV_EXIT_CODE 2
#define XO_EXIT_CODE  3

#endif
