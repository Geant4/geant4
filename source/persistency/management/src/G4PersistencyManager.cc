// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistencyManager.cc,v 2.7 1998/11/09 17:22:07 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4PersistencyManager 
//
// Implementation for concrete G4PersistencyManager.
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.06.20 Y.Morita  Implement geometry Store and Retrieve

#include "G4PersistencyManager.hh"

#include "G4ios.hh"

G4PersistencyManager* G4PersistencyManager::f_PersistencyManager = NULL;

G4PersistencyManager* G4PersistencyManager::GetPersistencyManager()
{
  if(!f_PersistencyManager)
  {
    f_PersistencyManager = new G4PersistencyManager;
  }
  return f_PersistencyManager;
}

G4PersistencyManager* G4PersistencyManager::get_PersistencyManagerIfExist()
{ return f_PersistencyManager; }

G4PersistencyManager::G4PersistencyManager()
{

  dbApp = new HepDbApplication("g4example");

  if (!dbApp)
  {
    G4cerr << "could not allocate HepDbApplication in G4PersistencyManager!" << endl;
  }

  // Open and Initialize a database

  // Note: The name of the federated database file "G4EXAMPLE" is 
  //       hard-coded here for now.
  //       In near future G4PersistencyManager should take an argument
  //       of the file name so that user can specify it on the fly.

  const G4String bootFileName = "G4EXAMPLE";
  G4cout << "Opening Federated Database " << bootFileName << "." << endl;
  dbApp->fdBootName(bootFileName);
  dbApp->Init();

  // Locate or create run databases and containers
  const G4String runDBName = "Runs";
  if( ! f_pRunMan.LocateDB( dbApp, runDBName) )
  {
    G4cerr << "could not locate run database in G4PersistencyManager!" << endl;
  }

  // Locate or create event databases and containers
  const G4String eventDBName = "Events";
  if( ! f_pEventMan.LocateDB( dbApp, eventDBName) )
  {
    G4cerr << "could not locate event database in G4PersistencyManager!" << endl;
  }

}

G4PersistencyManager::~G4PersistencyManager()
{
  // Any persistent objects which are created during save() should
  // not be deleted in ~G4PersistencyManager() implicitly, because
  // they are "persistent"!

  // On the contrary, any transient objects which are created during
  // Retrieve() should be deleted explicitly when they are no longer
  // in use.  It is the responsibility of the caller of Retrieve()
  // to destroy the retrieved transient objects if Retrieve() is to be
  // called repeatedlly.

  G4cout << "PersistencyManager is deleting." << endl;
}

//----------------------------------------------------------------------------

G4bool G4PersistencyManager::Store(const G4Event* anEvent)
{
  return f_pEventMan.Store(dbApp, anEvent);
}

G4bool G4PersistencyManager::Store(const G4Run* aRun)
{
  return f_pRunMan.Store(dbApp, aRun);
}

G4bool G4PersistencyManager::Store(const G4VPhysicalVolume* aWorld)
{
  return f_pGeomMan.Store(dbApp, aWorld);
}

G4bool G4PersistencyManager::Retrieve(G4Event*& anEvent)
{
  return f_pEventMan.Retrieve(dbApp, anEvent);
}

G4bool G4PersistencyManager::Retrieve(G4Run*& aRun)
{
  return f_pRunMan.Retrieve(dbApp, aRun);
}

G4bool G4PersistencyManager::Retrieve(G4VPhysicalVolume*& theWorld)
{
  return f_pGeomMan.Retrieve(dbApp, theWorld);
}

