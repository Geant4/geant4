// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentEventMan.hh,v 1.1 1999/01/07 16:10:56 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PersistentEventMan 
//
// A Utility class for storing and retrieving the event objects.
//
// This class is not persistent-capable. 
//
// Member functions:
// =================
//  G4bool Store(const G4Event* anEvent);
//    Store anEvent and associated objects into database.
//  G4bool Retrieve(G4Event*& anEvent);
//    Retrieve anEvent and associated objects into database.
//  G4Bool LocateDB( HepDBApplication* dbApp,
//                   const G4String eventDBName );
//    Create and locate event database and container.
//
// Member data:
// ============
//  HepDatabaseRef f_EventDB;
//    Reference to the event database.
//  HepContainerRef f_EventContainer;
//    Reference to the event container in the database.
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4PERSISTENTEVENTMAN_H
#define G4PERSISTENTEVENTMAN_H 1

#include "globals.hh"

#include "HepODBMS/clustering/HepDbApplication.h"

class G4Event;

class G4PersistentEventMan 
{
  public:
      G4PersistentEventMan();
      ~G4PersistentEventMan();
      HepDatabaseRef  GetDB();
      HepContainerRef GetContainer();
      G4bool LocateDB( HepDbApplication* dbApp,
                       const G4String eventDBName );

  private:
      HepDatabaseRef f_EventDB;
      HepContainerRef f_EventContainer;

  public:
      G4bool Store( HepDbApplication* dbApp,
                    const G4Event* anEvent );
      G4bool Retrieve( HepDbApplication* dbApp,
                       G4Event*& anEvent );

};


#endif

