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
// class G4AssemblyStore
//
// Class description:
//
// Container for all assemblies, with functionality derived from
// std::vector<T>. The class is a `singleton', in that only
// one can exist, and access is provided via the static method
// G4AssemblyStore::GetInstance().
//
// All assemblies should be registered with G4AssemblyStore, and removed on
// their destruction. The underlying container initially has a capacity of 20.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>
//
// Member data:
//
// static G4AssemblyStore*
//   - Pointer to the single G4AssemblyStore

// History:
// 9.10.18 G.Cosmo Initial version
// --------------------------------------------------------------------
#ifndef G4ASSEMBLYSTORE_HH
#define G4ASSEMBLYSTORE_HH

#include <vector>
#include "G4Types.hh"
#include "G4String.hh"
#include "G4VStoreNotifier.hh"

class G4AssemblyVolume;

class G4AssemblyStore : public std::vector<G4AssemblyVolume*>
{
  public:  // with description

    static void Register(G4AssemblyVolume* pAssembly);
      // Add the assembly to the collection.
    static void DeRegister(G4AssemblyVolume* pAssembly);
      // Remove the assembly from the collection.
    static G4AssemblyStore* GetInstance();
      // Get a ptr to the unique G4AssemblyStore, creating it if necessary.
    static void SetNotifier(G4VStoreNotifier* pNotifier);
      // Assign a notifier for allocation/deallocation of assemblies.
    static void Clean();
      // Delete all assemblies from the store.

    G4AssemblyVolume* GetAssembly(unsigned int Id, G4bool verbose=true) const;
      // Returns an assembly through its Id number specification.

  protected:

    G4AssemblyStore();
      // Protected singleton constructor.
    virtual ~G4AssemblyStore();
      // Destructor: takes care to delete allocated assemblies.

  private:

    static G4AssemblyStore* fgInstance;
    static G4ThreadLocal G4VStoreNotifier* fgNotifier;
    static G4ThreadLocal G4bool locked;
};

#endif
