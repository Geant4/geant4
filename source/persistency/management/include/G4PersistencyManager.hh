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
// $Id: G4PersistencyManager.hh,v 1.16 2001/07/11 10:02:26 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.06.20 Y.Morita  Implement geometry Store and Retrieve
// 98.10.30 Y.Morita  Splitted into event/run/geometry utility classes
// 99.10.22 Y.Morita  Allow constructor with the name of Boot file
// 99.11.12 Y.Morita  Splitted classes into separate objects

#ifndef G4PERSISTENCYMANAGER_HH
#define G4PERSISTENCYMANAGER_HH 1

#include "globals.hh"
#include "G4VPersistencyManager.hh"
#include "G4DatabaseTypes.hh"

#include "HepODBMS/odbms/HepODBMS.h"

// forward declarations
class G4PersistentHitMan;
class G4PersistentDigitMan;
class G4PersistentEventMan;
class G4PersistentRunMan;
class G4PersistentGeomMan;
class G4PersistencyMessenger;
class G4TransactionManager;
class HepDbApplication;

// Class Description:
//   A Class responsible for storing and retrieving the run, event,
//   Hit and geometry objects into ODBMS using HepODBMS interface.
//
//   The class is `singleton', with access via
//     G4PersistencyManager::GetPersistencyManager().
//
//   This class itself is not persistent-capable. 

class G4PersistencyManager 
 : public G4VPersistencyManager
{
  friend class G4PersistencyMessenger;

  public: // With description
      static G4PersistencyManager* GetPersistencyManager();
        // Return pointer to singleton instance of the class.
        // An instance is created if it does not exist.
      static G4PersistencyManager* get_PersistencyManagerIfExist();
        // Return pointer to singleton instance of the class.

  protected:
      // use GetPersistencyManager() instead
      G4PersistencyManager();

  public:
      ~G4PersistencyManager();

  private:
      static G4PersistencyManager* f_PersistencyManager;
      G4int  f_verboseLevel;

  public: // With description
      G4bool Store(const G4Event* anEvent);
        // stores anEvent and the associated objects into database.
      G4bool Store(const G4Run* aRun);
        // stores aRun and the associated objects into database.
      G4bool Store(const G4VPhysicalVolume* aWorld);
        // stores the world volume and the entire geometry informaion.

      G4bool Retrieve(G4Event*& anEvent);
        // retrieves anEvent and the associated objects from database.
      G4bool Retrieve(G4Run*& aRun);
        // retrieves aRun and the associated objects from database.
      G4bool Retrieve(G4VPhysicalVolume*& aWorld);
        // retrieves the world volume and the entire geometry informaion.

  protected:
      // protected interface (delegated to TransactionManager)
      HepDbApplication* DbApp();

  protected:
      // interface with PersistencyMessenger (delegated to TransactionManager)
      G4bool SelectDB(ETypeOfDB dbtype, G4String dbname, G4bool updateMode);

      G4String DBName(ETypeOfDB dbtype);
      G4String ContainerName(ETypeOfDB dbtype);
      HepDatabaseRef DBref(ETypeOfDB dbtype);
      HepContainerRef ContainerRef(ETypeOfDB dbtype);

      G4String DBContainerName(ETypeOfDB dbtype);

  public: // With description
      G4bool StartTransaction( ETypeOfDB dbtype,
                               ETransactionMode dbmode,
                               G4bool isSustained);
        // start an database transaction.
        //   ETypeOfDB:  kEventDB, kRunDB, kGeomDB
        //   ETransactionMode:  kUpdate, kRead
        //   isSustained = true  : transaction is sustained
        //                 false : transaction is atomic
        // See G4DatabaseTypes.hh for the definition.
      G4bool Commit(ETypeOfDB dbtype, G4bool isSustained);
        // commit the database transaction.
      G4bool  Abort(ETypeOfDB dbtype, G4bool isSustained);
        // abort the database transaction.

  private:
      // interface with PersistencyMessenger
      void SetVerboseLevel(G4int vl);
      inline G4int GetVerboseLevel()
      { return f_verboseLevel; };
      inline G4TransactionManager* GetTransactionManager()
      { return f_transactionMan; };

// Utility Classes
  private:
      G4PersistentHitMan*     f_pHitMan;
      G4PersistentDigitMan*   f_pDigitMan;
      G4PersistentEventMan*   f_pEventMan;
      G4PersistentRunMan*     f_pRunMan;
      G4PersistentGeomMan*    f_pGeomMan;
      G4PersistencyMessenger* f_persMessenger;
      G4TransactionManager*   f_transactionMan;
};

#endif

