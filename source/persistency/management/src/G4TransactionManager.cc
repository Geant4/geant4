//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TransactionManager.cc,v 1.8 2001/07/11 10:02:26 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// class G4TransactionManager 
//
// Implementation for concrete G4TransactionManager.
//
// History:
// 99.11.25 Y.Morita  Initial version

#include "G4TransactionManager.hh"

#include "G4ios.hh"

// forward declarations
#include "G4PersistentRunMan.hh"
#include "G4PersistentEventMan.hh"
#include "G4PersistentHitMan.hh"
#include "G4PersistentDigitMan.hh"
#include "G4PersistentGeomMan.hh"

G4TransactionManager::G4TransactionManager(
                           G4PersistentRunMan*   runMan,
                           G4PersistentEventMan* eventMan,
                           G4PersistentHitMan*   hitMan,
                           G4PersistentDigitMan* digitMan,
                           G4PersistentGeomMan*  geomMan)
 : f_RunDBName("Runs"), f_EventDBName("Events"), f_HitDBName("Events"),
   f_DigitDBName("Events"), f_GeomDBName("Geometry"),
   f_RunContainerName("RunContainer"), f_EventContainerName("EventContainer"),
   f_HitContainerName("EventContainer"), f_DigitContainerName("EventContainer"),
   f_GeomContainerName("GeomContainer"),
   f_pRunMan(runMan), f_pEventMan(eventMan), f_pHitMan(hitMan),
   f_pDigitMan(digitMan), f_pGeomMan(geomMan)
{
  const G4String applicationName = "g4example";
  f_dbApp = new HepDbApplication(applicationName);

  // Open and Initialize a database
  G4cout << "Opening federated database OO_FD_BOOT." << G4endl;
  f_dbApp->init();

  // Locate or create run,event,geometry databases and containers
  SelectDB( kRunDB,   f_RunDBName,   kUpdate);
  SelectDB( kEventDB, f_EventDBName, kUpdate);
  SelectDB( kGeomDB,  f_GeomDBName,  kUpdate);
}

G4TransactionManager::~G4TransactionManager()
{
  delete f_dbApp;
}

//----------------------------------------------------------------------------

G4bool G4TransactionManager::StartTransaction(
                                  ETypeOfDB dbtype,
                                  ETransactionMode dbmode,
                                  G4bool isSustained)
{
  ESustainedState f_sustainedState = CheckState(dbtype, isSustained);

  switch(f_sustainedState)
  {
    case kAlreadySelected:
      if(f_verboseLevel>2)
        G4cout << "G4TransactionManager: Sustained transaction for /"
               << DBName(f_whichDB) << "/ already in progress." << G4endl;
      break;
    case kCannotOverride:
      G4cerr << "G4TransactionManager: Cannot override the existing "
             << "transaction for /" << DBName(f_whichDB) << "/" << G4endl;
      return false;
      break;
    case kStartNewSustained:
      f_isSustained = isSustained;
      f_whichDB = dbtype;
      f_transactionMode = dbmode;
      DoStart(dbtype, dbmode, true);
      break;
    case kCommitAndStartNonSustained:
      f_dbApp->commit();
      if(f_verboseLevel>1)
        G4cout << "Sustained transaction for /" << DBName(f_whichDB)
               << "/ is paused." << G4endl;
      DoStart(dbtype, dbmode, false);
      break;
    case kStartNonSustained:
      DoStart(dbtype, dbmode, false);
      break;
  }

  return true;
}

ESustainedState G4TransactionManager::CheckState(
                                  ETypeOfDB dbtype, G4bool isSustained)
{
  if( isSustained )
    if( f_isSustained )
      if( f_whichDB == dbtype )
        return kAlreadySelected;
      else
        return kCannotOverride;
    else
      return kStartNewSustained;
  else
    if( f_isSustained )
      if( f_whichDB == dbtype )
        return kAlreadySelected;
      else
        return kCommitAndStartNonSustained;
    else
      return kStartNonSustained;
}

G4bool G4TransactionManager::DoStart(ETypeOfDB dbtype,
                                     ETransactionMode dbmode,
                                     G4bool isSustained)
{
  switch(dbmode)
  {
    case kUpdate:
      f_dbApp->startUpdate();
      break;
    case kRead:
      f_dbApp->startRead();
      break;
  }

  HepDatabaseRef aDBref = f_dbApp->db( DBName(dbtype) );
  if( aDBref == 0 )
  {
    G4cerr << "G4TransactionManager: Could not create or find /"
           << DBName(dbtype) << "/ database." << G4endl;
    return false;
  }
  SetDB( dbtype, aDBref );

  HepContainerRef aContRef = f_dbApp->container(ContainerName(dbtype));
  if( aContRef == 0 )
  {
    G4cerr << "G4TransactionManager: Could not create or find /"
           << ContainerName(dbtype) << "/ container." << G4endl;
    return false;
  }
  SetContainer( dbtype, aContRef );

  f_currentDB = dbtype;

  if(f_verboseLevel>1)
  {
    G4cout << "Transaction started for /" << DBName(dbtype) << "/"
           << ContainerName(dbtype) << "/ with ";
    if(isSustained)
      G4cout << "sustained mode." << G4endl;
    else
      G4cout << "non-sustained mode." << G4endl;
  }

  return true;
}

G4bool G4TransactionManager::Commit(ETypeOfDB dbtype, G4bool isSustained)
{
  if( dbtype == f_currentDB )
  {
    if( ! f_isSustained || isSustained )
    {
      f_dbApp->commit();

      if(f_verboseLevel>1)
        G4cout << "Transaction is committed on /" << DBName(dbtype)
               << "/" << ContainerName(dbtype) << "/" << G4endl;

      if( f_isSustained && dbtype != f_whichDB )
      {
        StartTransaction( f_whichDB, f_transactionMode, true );
        if(f_verboseLevel>1)
          G4cout << "Resumeing transaction on /" << DBName(f_whichDB)
                 << "/" << ContainerName(f_whichDB) << "/" << G4endl;
      }
      else
      {
        f_isSustained = false;
      }

      return true;
    }
  }
  else
  {
    G4cerr << "Error: Commit() received on /" << DBName(dbtype)
           << "/" << ContainerName(dbtype) << "/" << G4endl
           << "       Current transaction is on /" << DBName(f_currentDB)
           << "/" << ContainerName(f_currentDB) << "/" << G4endl;
    return false;
  }
}

G4bool G4TransactionManager::Abort(ETypeOfDB dbtype, G4bool isSustained)
{
  if( dbtype == f_currentDB )
  {
    if( ! f_isSustained || isSustained )
    {
      f_dbApp->abort();

      if(f_verboseLevel>1)
        G4cout << "Transaction is aborted on /" << DBName(dbtype)
               << "/" << ContainerName(dbtype) << "/" << G4endl;

      if( f_isSustained && dbtype != f_whichDB )
      {
        StartTransaction( f_whichDB, f_transactionMode, true );
        if(f_verboseLevel>1)
          G4cout << "Resumeing transaction on /" << DBName(f_whichDB)
                 << "/" << ContainerName(f_whichDB) << "/" << G4endl;
      }
      else
      {
        f_isSustained = false;
      }

      return true;
    }
  }
  else
  {
    G4cerr << "Error: Abort() received on /" << DBName(dbtype)
           << "/" << ContainerName(dbtype) << "/" << G4endl
           << "       Current transaction is on /" << DBName(f_currentDB)
           << "/" << ContainerName(f_currentDB) << "/" << G4endl;
    return false;
  }
}

G4bool G4TransactionManager::SelectDB( ETypeOfDB dbtype,
                                       G4String dbname,
                                       G4bool updateMode)
{
  G4bool theStatus = false;

  G4String oldDBname = DBName(dbtype);
  SetDBName(dbtype, dbname);

  // start a non-sustained transaction for this database type
  if( StartTransaction( dbtype, kUpdate, false ) )
  {
    Commit( dbtype, false );
    // update database ref and container ref in SubDbManagers
    SendDB(dbtype);
    SendContainer(dbtype);
    theStatus = true;

    if(f_verboseLevel>0)
      G4cout << "Set database to /"
             << DBName(dbtype) << "/." << G4endl;
  }
  else
  {
    Abort( dbtype, false );
    SetDBName(dbtype, oldDBname);
    G4cerr << "G4TransactionManager: Failed to set database /"
             << DBName(dbtype) << "/." << G4endl;
    theStatus = false;
  }

  return theStatus;
}

G4String G4TransactionManager::DBName(ETypeOfDB dbtype)
{
  switch(dbtype)
  {
    case kRunDB:
      return f_RunDBName;
      break;
    case kEventDB:
      return f_EventDBName;
      break;
    case kGeomDB:
      return f_GeomDBName;
      break;
  }
}

G4String G4TransactionManager::ContainerName(ETypeOfDB dbtype)
{
  switch(dbtype)
  {
    case kRunDB:
      return f_RunContainerName;
      break;
    case kEventDB:
      return f_EventContainerName;
      break;
    case kGeomDB:
      return f_GeomContainerName;
      break;
  }
}

HepDatabaseRef G4TransactionManager::DBref(ETypeOfDB dbtype)
{
  switch(dbtype)
  {
    case kRunDB:
      return f_RunDB;
      break;
    case kEventDB:
      return f_EventDB;
      break;
    case kGeomDB:
      return f_GeomDB;
      break;
  }
}

HepContainerRef G4TransactionManager::ContainerRef(ETypeOfDB dbtype)
{
  switch(dbtype)
  {
    case kRunDB:
      return f_RunContainer;
      break;
    case kEventDB:
      return f_EventContainer;
      break;
    case kGeomDB:
      return f_GeomContainer;
      break;
  }
}

void G4TransactionManager::SetDBName(ETypeOfDB dbtype, G4String dbname)
{
  switch(dbtype)
  {
    case kRunDB:
      f_RunDBName = dbname;
      break;
    case kEventDB:
      f_EventDBName = dbname;
      f_HitDBName   = dbname;
      f_DigitDBName = dbname;
      break;
    case kGeomDB:
      f_GeomDBName = dbname;
      break;
  }
}

void G4TransactionManager::SetDB(ETypeOfDB dbtype, HepDatabaseRef aDB)
{
  switch(dbtype)
  {
    case kRunDB:
      f_RunDB = aDB;
      break;
    case kEventDB:
      f_EventDB = aDB;
      f_HitDB   = aDB;
      f_DigitDB = aDB;
      break;
    case kGeomDB:
      f_GeomDB = aDB;
      break;
  }
}

void G4TransactionManager::SetContainer(ETypeOfDB dbtype, HepContainerRef aCont)
{
  switch(dbtype)
  {
    case kRunDB:
      f_RunContainer = aCont;
      break;
    case kEventDB:
      f_EventContainer = aCont;
      f_HitContainer   = aCont;
      f_DigitContainer = aCont;
      break;
    case kGeomDB:
      f_GeomContainer = aCont;
      break;
  }
}

void G4TransactionManager::SendDB(ETypeOfDB dbtype)
{
  switch(dbtype)
  {
    case kRunDB:
      f_pRunMan->SetDB(f_RunDB);
      break;
    case kEventDB:
      f_pEventMan->SetDB(f_EventDB);
      f_pHitMan  ->SetDB(f_EventDB);
      f_pDigitMan->SetDB(f_EventDB);
      break;
    case kGeomDB:
      f_pGeomMan->SetDB(f_GeomDB);
      break;
  }
}

void G4TransactionManager::SendContainer(ETypeOfDB dbtype)
{
  switch(dbtype)
  {
    case kRunDB:
      f_pRunMan->SetContainer(f_RunContainer);
      break;
    case kEventDB:
      f_pEventMan->SetContainer(f_EventContainer);
      f_pHitMan  ->SetContainer(f_EventContainer);
      f_pDigitMan->SetContainer(f_EventContainer);
      break;
    case kGeomDB:
      f_pGeomMan->SetContainer(f_GeomContainer);
      break;
  }
}

G4String G4TransactionManager::DBContainerName(ETypeOfDB dbtype)
{
  G4String aString = "/" + DBName(dbtype) + "/" + ContainerName(dbtype) + "/";
  return aString;
}

