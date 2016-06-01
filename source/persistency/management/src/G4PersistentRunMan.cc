// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentRunMan.cc,v 2.2 1998/11/10 04:05:52 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4PersistentRunMan 
//
// Implementation for concrete G4PersistentRunMan.
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.06.20 Y.Morita  Implement geometry Store and Retrieve

#include "G4PersistentRunMan.hh"

#include "G4Run.hh"
#include "G4PRun.hh"

#include "G4ios.hh"

G4PersistentRunMan::G4PersistentRunMan()
{;}

G4PersistentRunMan::~G4PersistentRunMan()
{;}

HepDatabaseRef G4PersistentRunMan::GetDB()
{
  return f_RunDB;
}

HepContainerRef G4PersistentRunMan::GetContainer()
{
  return f_RunContainer;
}

G4PersistentRunMan::LocateDB(HepDbApplication* dbApp,
                               const G4String runDBName)
{
  G4bool theStatus = true;

  // start a transaction to create databases and containers
  dbApp->startUpdate();

  // create a root named object tree namespace
//  namedRoot = new NamedNode;
//  NamedObjectTree  namedTree( namedRoot );

  // use database RunDB
  f_RunDB = dbApp->db("RunDB");

  if (f_RunDB == NULL)
  {
    G4cerr << "Could not create or find RunDB database." << endl;
    theStatus = false;
  }

  // create a new container for Run collection in this database
  f_RunContainer = dbApp->container("RunContainer"); 

  if ( f_RunContainer == NULL )
  {
    G4cerr << "could not find or create RunContainer in the database" << endl;
    theStatus = false;
  }

  // create a "RunDB" named object tree in the namespace
//  namedTree.mkdir("RunDB");
//  namedTree.cd("RunDB");

  // Commit the update
  dbApp->commit();

  return theStatus;
}

//----------------------------------------------------------------------------

G4bool G4PersistentRunMan::Store( HepDbApplication* dbApp,
                                  const G4Run* aRun)
{
// Start Update G4PRun in the database
  dbApp->startUpdate();

// Create persistent Run object
  new(f_RunContainer) G4PRun(aRun);

// Commit the update
  dbApp->commit();

#ifndef G4PERSISTENCY_DEBUG
  G4cout << "G4PRun " << aRun->GetRunID() << " Committed." << endl;
#endif

  return true;
}

G4bool G4PersistentRunMan::Retrieve( HepDbApplication* dbApp,
                                     G4Run*& aRun)
{
// not implemented yet

  G4bool fStatus = false;

  fStatus = true;
  return fStatus;
}

