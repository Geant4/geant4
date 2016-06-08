// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentRunMan.hh,v 1.1 1999/01/07 16:10:56 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PersistentRunMan 
//
// A Utility class for storing and retrieving the run objects.
//
// This class is not persistent-capable. 
//
// Member functions:
// =================
//  G4bool Store(const G4Run* aRun);
//    Store aRun and associated objects into database.
//  G4bool Retrieve(G4Run*& aRun);
//    Retrieve aRun and associated objects into database.
//  G4Bool LocateDB( HepDbApplication* dbApp,
//                   const G4String runDBName );
//    Create and locate run database and container.
//
// Member data:
// ============
//  HepDatabaseRef f_RunDB;
//    Reference to the run database.
//  HepContainerRef f_RunContainer;
//    Reference to the run container in the database.
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4PERSISTENTRUNMAN_H
#define G4PERSISTENTRUNMAN_H 1

#include "globals.hh"

#include "HepODBMS/clustering/HepDbApplication.h"

class G4Run;


class G4PersistentRunMan 
{
  public:
      G4PersistentRunMan();
      ~G4PersistentRunMan();
      HepDatabaseRef  GetDB();
      HepContainerRef GetContainer();
      G4bool LocateDB( HepDbApplication* dbApp,
                       const G4String runDBName );

  private:
      HepDatabaseRef f_RunDB;
      HepContainerRef f_RunContainer;

  public:
      G4bool Store( HepDbApplication* dbApp, const G4Run* aRun);
      G4bool Retrieve( HepDbApplication* dbApp, G4Run*& aRun );
};


#endif

