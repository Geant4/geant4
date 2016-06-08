// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransactionManager.hh,v 1.5.2.1 1999/12/07 20:50:14 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

// class description:
//
//      This is a class for handling database transaction
//      for G4PersistencyManager.
//
// Member functions:
// =================
//  G4TransactionManager();
//  ~G4TransactionManager();
//
// History:
// 99.11.25 Y.Morita  Initial version

#ifndef G4TransactionManager_HH
#define G4TransactionManager_HH 1

#include "globals.hh"
#include "G4VPersistentSubMan.hh"

// forward declarations
class G4PersistentRunMan;
class G4PersistentEventMan;
class G4PersistentHitMan;
class G4PersistentDigitMan;
class G4PersistentGeomMan;

#include "HepODBMS/clustering/HepDbApplication.h"

#include "G4DatabaseTypes.hh"

class G4TransactionManager 
 : public G4VPersistentSubMan
{
  friend class G4PersistencyManager;

  protected:
      G4TransactionManager(G4PersistentRunMan*   runMan,
                           G4PersistentEventMan* eventMan,
                           G4PersistentHitMan*   hitMan,
                           G4PersistentDigitMan* digitMan,
                           G4PersistentGeomMan*  geomMan);
      ~G4TransactionManager();

  private:
      HepDbApplication* f_dbApp;

      G4bool           f_isSustained;
      ETypeOfDB        f_whichDB;
      ETransactionMode f_transactionMode;
      ETypeOfDB        f_currentDB;

      HepDatabaseRef f_RunDB;
      HepDatabaseRef f_EventDB;
      HepDatabaseRef f_HitDB;
      HepDatabaseRef f_DigitDB;
      HepDatabaseRef f_GeomDB;

      G4String f_RunDBName;
      G4String f_EventDBName;
      G4String f_HitDBName;
      G4String f_DigitDBName;
      G4String f_GeomDBName;

      HepContainerRef f_RunContainer;
      HepContainerRef f_EventContainer;
      HepContainerRef f_HitContainer;
      HepContainerRef f_DigitContainer;
      HepContainerRef f_GeomContainer;

      const G4String f_RunContainerName;
      const G4String f_EventContainerName;
      const G4String f_HitContainerName;
      const G4String f_DigitContainerName;
      const G4String f_GeomContainerName;

  private:
      // interface with PersistencyManager
      inline HepDbApplication* DbApp()
      { return f_dbApp; }
      G4bool StartTransaction( ETypeOfDB dbtype,
                               ETransactionMode dbmode,
                               G4bool isSustained);
      G4bool Commit(ETypeOfDB dbtype, G4bool isSustained);
      G4bool  Abort(ETypeOfDB dbtype, G4bool isSustained);
      inline G4bool SustainedMode()
      { return f_isSustained; }

  private:
      // interface with PersistencyManager
      G4bool SelectDB(ETypeOfDB dbtype, G4String dbname, G4bool updateMode);

      G4String DBName(ETypeOfDB dbtype);
      G4String ContainerName(ETypeOfDB dbtype);
      HepDatabaseRef DBref(ETypeOfDB dbtype);
      HepContainerRef ContainerRef(ETypeOfDB dbtype);

      G4String DBContainerName(ETypeOfDB dbtype);

  private:
      // for internal use
      ESustainedState CheckState(ETypeOfDB dbtype, G4bool isSustained);
      G4bool DoStart(ETypeOfDB dbtype,
                     ETransactionMode dbmode,
                     G4bool isSustained);
      void SetDBName(ETypeOfDB dbtype, G4String dbname);
      void SetDB(ETypeOfDB dbtype, HepDatabaseRef aDB);
      void SetContainer(ETypeOfDB dbtype, HepContainerRef aCont);
      void SendDB(ETypeOfDB dbtype);
      void SendContainer(ETypeOfDB dbtype);

// Utility Classes
  private:
      G4PersistentRunMan*   f_pRunMan;
      G4PersistentEventMan* f_pEventMan;
      G4PersistentHitMan*   f_pHitMan;
      G4PersistentDigitMan* f_pDigitMan;
      G4PersistentGeomMan*  f_pGeomMan;
};


#endif

