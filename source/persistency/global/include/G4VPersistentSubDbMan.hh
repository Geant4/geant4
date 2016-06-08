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
// $Id: G4VPersistentSubDbMan.hh,v 1.5 2001/07/11 10:02:25 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $

// Class Description:
//   Abstract class to be used by G4PersistencyManager
//   and G4TransactionManager
//
// Member functions:
// =================
//  void SetDB(HepDatabaseRef aDB);
//    Set DB reference.
//  HepDatabaseRef GetDB();
//    Get DB reference.
//  void SetContainer(HepContainerRef container);
//    Set Container reference.
//  HepContainerRef GetContainer();
//    Get Container reference.
//
//  Note: inherited concrete class should implement Store() and Retrieve()
//        and also should set f_currentContainer in Retieve()
//    Store(HepDbApplication* dbApp, const object_type* object);
//    Retrieve(HepDbApplication* dbApp, object_type*& object);
//
// Member data:
// ============
//  HepContainerRef f_DB;
//    Reference to the database.
//  HepContainerRef f_container;
//    Reference to the container in the database.
//  HepContainerRef f_currentContainer;
//    Reference to the current scan scope in Retrieve().
//
//    example in a piece of Retrive():
//      if( f_container != f_currentContainer )
//      {
//        f_currentContainer = f_container;
//        pRun_iterator.scan(f_currentContainer);
//      }

// History:
// 99.11.25 Y.Morita  Initial creation

#ifndef G4VPersistentSubDbMan_hh
#define G4VPersistentSubDbMan_hh 1

#include "globals.hh"
#include "HepODBMS/odbms/HepODBMS.h"

#include "G4PersistentTypes.hh"
#include "G4DatabaseTypes.hh"

class G4VPersistentSubDbMan 
{
  friend class G4PersistencyManager;
  friend class G4TransactionManager;

  protected:
      // to be constructed by G4PersistencyManager only
      G4VPersistentSubDbMan();
      virtual ~G4VPersistentSubDbMan();

  protected:
      // interface with G4PersistencyManager
      // Store(HepDbApplication* dbApp, const object_type* object);
      // Retrieve(HepDbApplication* dbApp, object_type*& object);

  protected:
      // interface with G4TransactionManager
      inline void SetDB(HepDatabaseRef aDB)
      { f_DB = aDB; }
      inline void SetContainer(HepContainerRef container)
      { f_container = container; }

  public: // With description
      inline HepDatabaseRef GetDB()
      { return f_DB; }
        // Returns the reference pointer of the database
      inline HepContainerRef GetContainer()
      { return f_container; }
        // Returns the reference pointer of the container

  protected:
      HepDatabaseRef  f_DB;
      HepContainerRef f_container;
      HepContainerRef f_currentContainer;
};

#endif

