// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef G4VInteractiveSession_H
#define G4VInteractiveSession_H 1

#include <rw/tvhdict.h>

#include "G4VInteractorManager.hh"

class G4UImessenger;

class G4VInteractiveSession {
public:
  G4VInteractiveSession();
  virtual ~G4VInteractiveSession();
  virtual void AddMenu (const char*,const char*);
  virtual void AddButton (const char*,const char*,const char*);
  void AddInteractor(G4String,G4Interactor);
  G4Interactor GetInteractor(G4String);
private:
  G4UImessenger* messenger;
  RWTValHashDictionary<G4String,G4Interactor> interactors;
};

#endif

