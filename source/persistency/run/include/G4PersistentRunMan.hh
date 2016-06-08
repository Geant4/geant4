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
// $Id: G4PersistentRunMan.hh,v 1.10 2001/07/11 10:02:26 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

// Class Description:
//   A Utility class to be used by G4PersistencyManager.
// Average users do not need to use this class.
//
// This class is not persistent-capable.
//

// Member functions:
// =================
//  G4bool Store(HepDbApplication* dbApp, const G4Run* aRun);
//    Store aRun and associated objects into database.
//  G4bool Retrieve(HepDbApplication* dbApp, G4Run*& aRun);
//    Retrieve aRun and associated objects into database.
//  void SetContainer(HepContainerRef container);
//    Set Run Container reference.
//
// Member data:
// ============
//  HepContainerRef f_container;
//    Reference to the event container in the database.
//  HepContainerRef f_currentContainer;
//    Reference to the current scan scope in Retrieve().
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4PERSISTENTRUNMAN_H
#define G4PERSISTENTRUNMAN_H 1

#include "globals.hh"

#include "G4VPersistentSubMan.hh"
#include "G4VPersistentSubDbMan.hh"
#include "G4PRun.hh"

class HepDbApplication;
class G4Run;

class G4PersistentRunMan
 : public G4VPersistentSubMan, public G4VPersistentSubDbMan
{
  friend class G4PersistencyManager;
  friend class G4TransactionManager;

  private:
      G4PersistentRunMan();
      ~G4PersistentRunMan();

  private:
      G4int f_currentRunID;
      HepRef(G4PRun) f_currentPRun;

  public:
      inline G4int CurrentRunID()
      { return f_currentRunID; }

  protected:
      // interface with G4TransactionManager

  private:
      // interface with G4PersistencyManager
      G4bool Store( HepDbApplication* dbApp, const G4Run* aRun);
      G4bool Retrieve( HepDbApplication* dbApp, G4Run*& aRun );
};


#endif

