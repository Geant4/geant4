// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentEventMan.cc,v 2.1 1998/11/09 17:19:46 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
// class G4PersistentEventMan 
//
// Implementation for concrete G4PersistentEventMan.
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.06.20 Y.Morita  Implement geometry Store and Retrieve

#include "G4PersistentEventMan.hh"

#include "G4Event.hh"
#include "G4PEvent.hh"

#include "G4ios.hh"

G4PersistentEventMan::G4PersistentEventMan()
{;}

G4PersistentEventMan::~G4PersistentEventMan()
{;}

HepDatabaseRef G4PersistentEventMan::GetDB()
{
  return f_EventDB;
}

HepContainerRef G4PersistentEventMan::GetContainer()
{
  return f_EventContainer;
}

G4bool G4PersistentEventMan::LocateDB(HepDbApplication* dbApp,
                                      const G4String eventDBName)
{
  G4bool theStatus = true;

  // start a transaction to create databases and containers
  dbApp->startUpdate();

  // create a root named object tree namespace
//  namedRoot = new NamedNode;
//  NamedObjectTree  namedTree( namedRoot );

  // use database Events
  f_EventDB = dbApp->db(eventDBName);

  if (f_EventDB == NULL)
  {
    G4cerr << "Could not create or find " << eventDBName << " database." << endl;
    theStatus = false;
  }

  // create a new container for event collection in this database
  f_EventContainer = dbApp->container("EventContainer"); 

  if ( f_EventContainer == NULL )
  {
    G4cerr << "could not find or create eventContainer in the database" << endl;
    theStatus = false;
  }

  // create a "Events" named object tree in the namespace
//  namedTree.mkdir(eventDBName);
//  namedTree.cd(eventDBName);

  // Commit the update
  dbApp->commit();

  return theStatus;
}

//----------------------------------------------------------------------------

G4bool G4PersistentEventMan::Store( HepDbApplication* dbApp,
                                    const G4Event* anEvent )
{
// Start Update G4PEvent in the database
  dbApp->startUpdate();

// Create persistent event object
  new(f_EventContainer) G4PEvent(anEvent);

// Put it into event iterator

// Commit the update
  dbApp->commit();

#ifndef G4PERSISTENCY_DEBUG
  G4cout << "G4PEvent " << anEvent->GetEventID() << " Committed." << endl;
#endif

  return true;
}

G4bool G4PersistentEventMan::Retrieve( HepDbApplication* dbApp,
                                       G4Event*& anEvent )
{
// not implemented yet

  G4bool fStatus = false;

  // loop thru the event iterator

  fStatus = true;
  return fStatus;
}

