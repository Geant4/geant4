// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistencyManager.cc,v 1.16 1999/12/15 14:51:27 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
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

// forward declarations
#include "G4PersistentHitMan.hh"
#include "G4PersistentDigitMan.hh"
#include "G4PersistentEventMan.hh"
#include "G4PersistentRunMan.hh"
#include "G4PersistentGeomMan.hh"
#include "G4PersistencyMessenger.hh"
#include "G4TransactionManager.hh"
#include "HepODBMS/clustering/HepDbApplication.h"

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
 : f_verboseLevel(0)
{
  f_pHitMan   = G4PersistentHitMan::GetPersistentHitMan();
  f_pDigitMan = G4PersistentDigitMan::GetPersistentDigitMan();
  f_pEventMan = new G4PersistentEventMan(f_pHitMan, f_pDigitMan);
  f_pRunMan   = new G4PersistentRunMan;
  f_pGeomMan  = new G4PersistentGeomMan;

  f_persMessenger  = new G4PersistencyMessenger(this);
  f_transactionMan = new G4TransactionManager
                          (f_pRunMan, f_pEventMan, f_pHitMan, f_pDigitMan, 
                           f_pGeomMan);
}

G4PersistencyManager::~G4PersistencyManager()
{
  // delete utility classes

  delete f_pHitMan;
  delete f_pDigitMan;
  delete f_pEventMan;
  delete f_pRunMan;
  delete f_pGeomMan;
  delete f_persMessenger;
  delete f_transactionMan;

  // Any persistent objects which are created during Store()
  // should not be deleted in ~G4PersistencyManager(), because
  // they are "persistent"!

  if(f_verboseLevel>0)
    G4cout << "PersistencyManager is deleting." << G4endl;
}

//----------------------------------------------------------------------------

G4bool G4PersistencyManager::Store(const G4Event* anEvent)
{
  f_transactionMan->StartTransaction(kEventDB, kUpdate, false);

  if( f_pEventMan->Store(f_transactionMan->DbApp(), anEvent) )
  {
    if( f_verboseLevel>1 )
    {
      G4cout << " -- G4PEvent " << f_pEventMan->CurrentEventID() << " stored in "
             << DBContainerName(kEventDB);
      if( f_transactionMan->SustainedMode() )
        G4cout << " (sustained)" << G4endl;
      else
        G4cout << G4endl;
    }
    f_transactionMan->Commit(kEventDB, false);
    return true;
  }
  else
  {
    f_transactionMan->Abort(kEventDB, false);
    G4cerr << "G4PersistencyManager: Failed to store G4PEvent in "
           << DBContainerName(kEventDB) << G4endl;
    return false;
  }
}

G4bool G4PersistencyManager::Store(const G4Run* aRun)
{
  f_transactionMan->StartTransaction(kRunDB, kUpdate, false);

  if( f_pRunMan->Store(f_transactionMan->DbApp(), aRun) )
  {
    if( f_verboseLevel>0 )
      G4cout << " -- G4PRun " << f_pRunMan->CurrentRunID() << " stored in "
             << DBContainerName(kRunDB) << G4endl;
    f_transactionMan->Commit(kRunDB, false);
    return true;
  }
  else
  {
    f_transactionMan->Abort(kRunDB, false);
    G4cerr << "G4PersistencyManager: Failed to store G4PRun in "
           << DBContainerName(kRunDB) << G4endl;
    return false;
  }
}

G4bool G4PersistencyManager::Store(const G4VPhysicalVolume* aWorld)
{
  f_transactionMan->StartTransaction(kGeomDB, kUpdate, false);

  if( f_pGeomMan->Store(f_transactionMan->DbApp(), aWorld) )
  {
    if( f_verboseLevel>0 )
      G4cout << " -- Geometry stored in "
             << DBContainerName(kGeomDB) << G4endl;
    f_transactionMan->Commit(kGeomDB, false);
    return true;
  }
  else
  {
    f_transactionMan->Abort(kGeomDB, false);
    G4cerr << "G4PersistencyManager: Failed to store Geometry in "
           << DBContainerName(kGeomDB) << G4endl;
    return false;
  }
}

G4bool G4PersistencyManager::Retrieve(G4Event*& anEvent)
{
  f_transactionMan->StartTransaction(kEventDB, kRead, false);

  if( f_pEventMan->Retrieve(f_transactionMan->DbApp(), anEvent) )
  {
    if( f_verboseLevel>1 )
    {
      if( anEvent )
        G4cout << " -- G4Event " << f_pEventMan->CurrentEventID()
               << " retrieved from "
               << DBContainerName(kEventDB) << G4endl;
      else
        G4cout << " -- scan of G4Event from "
               << DBContainerName(kEventDB)
               << " is completed." << G4endl;
    }
    f_transactionMan->Commit(kEventDB, false);
    return true;
  }
  else
  {
    f_transactionMan->Abort(kEventDB, false);
    G4cerr << "G4PersistencyManager: Failed to retrieve G4Event from "
           << DBContainerName(kEventDB) << G4endl;
    return false;
  }
}

G4bool G4PersistencyManager::Retrieve(G4Run*& aRun)
{
  f_transactionMan->StartTransaction(kRunDB, kRead, false);

  if( f_pRunMan->Retrieve(f_transactionMan->DbApp(), aRun) )
  {
    if( f_verboseLevel>0 )
    {
      if( aRun )
        G4cout << " -- G4Run " << f_pRunMan->CurrentRunID()
               << " retrieved from "
               << DBContainerName(kRunDB) << "." << G4endl;
      else
        G4cout << " -- scan of G4Run from "
               << DBContainerName(kRunDB)
               << " is completed." << G4endl;
    }
    f_transactionMan->Commit(kRunDB, false);
    return true;
  }
  else
  {
    f_transactionMan->Abort(kRunDB, false);
    G4cerr << "G4PersistencyManager: Failed to retrieve G4Run from "
           << DBContainerName(kRunDB) << G4endl;
    return false;
  }
}

G4bool G4PersistencyManager::Retrieve(G4VPhysicalVolume*& theWorld)
{
  f_transactionMan->StartTransaction(kGeomDB, kRead, false);

  if( f_pGeomMan->Retrieve(f_transactionMan->DbApp(), theWorld) )
  {
    if( f_verboseLevel>0 )
    {
      if( theWorld )
        G4cout << " -- Geometry retrieved from "
               << DBContainerName(kGeomDB) << G4endl;
    }
    f_transactionMan->Commit(kGeomDB, false);
    return true;
  }
  else
  {
    f_transactionMan->Abort(kGeomDB, false);
    G4cerr << "G4PersistencyManager: Failed to retrieve Geometry from "
           << DBContainerName(kGeomDB) << G4endl;
    return false;
  }
}


//----------------------------------------------------------------------------

HepDbApplication* G4PersistencyManager::DbApp()
{ return f_transactionMan->DbApp(); }

//----------------------------------------------------------------------------

G4bool G4PersistencyManager::SelectDB(ETypeOfDB dbtype,
                                      G4String dbname, G4bool updateMode)
{ return f_transactionMan->SelectDB(dbtype, dbname, updateMode); }

G4String G4PersistencyManager::DBName(ETypeOfDB dbtype)
{ return f_transactionMan->DBName(dbtype); }

G4String G4PersistencyManager::ContainerName(ETypeOfDB dbtype)
{ return f_transactionMan->ContainerName(dbtype); }

HepDatabaseRef G4PersistencyManager::DBref(ETypeOfDB dbtype)
{ return f_transactionMan->DBref(dbtype); }

HepContainerRef G4PersistencyManager::ContainerRef(ETypeOfDB dbtype)
{ return f_transactionMan->ContainerRef(dbtype); }

G4String G4PersistencyManager::DBContainerName(ETypeOfDB dbtype)
{ return f_transactionMan->DBContainerName(dbtype); }

//----------------------------------------------------------------------------

G4bool G4PersistencyManager::StartTransaction( ETypeOfDB dbtype,
                                               ETransactionMode dbmode,
                                               G4bool isSustained)
{
  return f_transactionMan->StartTransaction( 
                                      dbtype, dbmode, isSustained );
}

G4bool G4PersistencyManager::Commit(ETypeOfDB dbtype, G4bool isSustained)
{ return f_transactionMan->Commit(dbtype, isSustained); }

G4bool G4PersistencyManager::Abort(ETypeOfDB dbtype, G4bool isSustained)
{ return f_transactionMan->Abort(dbtype, isSustained); }

//----------------------------------------------------------------------------

void G4PersistencyManager::SetVerboseLevel(G4int vl)
{
  f_verboseLevel = vl;
  f_pHitMan->SetVerboseLevel(vl);
  f_pDigitMan->SetVerboseLevel(vl);
  f_pRunMan->SetVerboseLevel(vl);
  f_pEventMan->SetVerboseLevel(vl);
  f_pGeomMan->SetVerboseLevel(vl);
  f_transactionMan->SetVerboseLevel(vl);
}

