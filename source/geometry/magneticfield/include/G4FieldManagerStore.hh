//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4FieldManagerStore.hh 103228 2017-03-22 14:52:32Z gcosmo $
//
// class G4FieldManagerStore
//
// Class description:
//
// Container for all FieldManagers, with functionality derived from
// std::vector<T>. The class is a `singleton', in that only
// one can exist, and access is provided via the static function
// G4FieldManagerStore::GetInstance()
//
// All FieldManagers should be registered with G4FieldManagerStore,
// and removed on their destruction. 
// Intended principally to enable resetting of 'state' at start of event.
// The underlying container initially has a capacity of 100.
//
// Member data:
//
// static G4FieldManagerStore* fgInstance
//   - Ptr to the single G4FieldManagerStore.

// History:
// 07.12.07 J.Apostolakis  Initial version
// --------------------------------------------------------------------
#ifndef G4FIELDMANAGERSTORE_HH
#define G4FIELDMANAGERSTORE_HH

#include <vector>

#include "G4FieldManager.hh"

class G4FieldManagerStore : public std::vector<G4FieldManager*>
{
  public:  // with description

    static void Register(G4FieldManager* pVolume);
      // Add the logical volume to the collection.
    static void DeRegister(G4FieldManager* pVolume);
      // Remove the logical volume from the collection.
    static G4FieldManagerStore* GetInstance();
      // Get a ptr to the unique G4FieldManagerStore, creating it if necessary.
    static G4FieldManagerStore* GetInstanceIfExist();
      // Get a ptr to the unique G4FieldManagerStore.
    static void Clean();
      // Delete all volumes from the store.

    void ClearAllChordFindersState();
      // Looping over all field managers, call each one to reset step estimate

    ~G4FieldManagerStore();
      // Destructor: takes care to delete allocated field managers.

  protected:
  
    G4FieldManagerStore();

  private:

    static G4ThreadLocal G4FieldManagerStore* fgInstance;
    static G4ThreadLocal G4bool locked;
};

#endif
