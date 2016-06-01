// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistencyManager.hh,v 2.5 1998/11/09 17:21:42 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4PersistencyManager 
//
// A Class responsible for storing and retrieving the run, event,
// Hit and geometry objects into ODBMS using HepODBMS interface.
//
// The class is `singleton', with access via
//   G4PersistencyManager::GetPersistencyManager.
//
// This class itself is not persistent-capable. 
//
// Member functions:
// =================
//  static G4PersistencyManager* GetPersistencyManager();
//  static G4PersistencyManager* get_PersistencyManagerIfExist();
//   Return ptr to singleton instance of the class.
//
//  G4bool Store(const G4Event* anEvent);
//    Store anEvent and associated objects into a database.
//  G4bool Store(const G4Run* aRun);
//    Store aRun and associated objects into a database.
//  G4bool Store(const G4VPhysicalVolume* aWorld);
//    Store aWorld and entire geometry object tree into a database.
//
//  G4bool Retrieve(G4Event*& anEvent);
//    Retrieve anEvent and associated objects from a database.
//  G4bool Retrieve(G4Run*& aRun);
//    Retrieve aRun and associated objects from a database.
//  G4bool Retrieve(G4VPhysicalVolume*& aWorld);
//    Retrieve aWorld and entire geometry object tree from a database.
//
// Member data:
// ============
// static G4PersistencyManager* fPersistencyManager
//   Ptr to the unique instance of class
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.06.20 Y.Morita  Implement geometry Store and Retrieve
// 98.10.30 Y.Morita  Splitted into event/run/geometry utility classes

#ifndef G4PERSISTENCYMANAGER_HH
#define G4PERSISTENCYMANAGER_HH 1

#include "globals.hh"
#include "G4VPersistencyManager.hh"

#include "G4PersistentEventMan.hh"
#include "G4PersistentRunMan.hh"
#include "G4PersistentGeomMan.hh"

////#include "HepODBMS/clustering/HepContainerHint.h"
#include "HepODBMS/clustering/HepDbApplication.h"

class G4PersistencyManager 
 : public G4VPersistencyManager
{
  public:
      static G4PersistencyManager* GetPersistencyManager();
      static G4PersistencyManager* get_PersistencyManagerIfExist();
      G4PersistencyManager();

  public:
      ~G4PersistencyManager();

  private:
      HepDbApplication* dbApp;

//      static HepRef(NamedNode) namedRoot;
//      static NamedObjectTree* 

  private: 
      static G4PersistencyManager * f_PersistencyManager;

  public:
      G4bool Store(const G4Event* anEvent);
      G4bool Store(const G4Run* aRun);
      G4bool Store(const G4VPhysicalVolume* aWorld);

      G4bool Retrieve(G4Event*& anEvent);
      G4bool Retrieve(G4Run*& aRun);
      G4bool Retrieve(G4VPhysicalVolume*& aWorld);

// Utility Classes
  private:
      G4PersistentEventMan f_pEventMan;
      G4PersistentRunMan f_pRunMan;
      G4PersistentGeomMan f_pGeomMan;

};


#endif

