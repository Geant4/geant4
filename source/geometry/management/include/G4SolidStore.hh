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
// $Id: G4SolidStore.hh,v 1.6 2001-07-11 09:59:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SolidStore
//
// Class description:
//
// Container for all solids, with functionality derived from
// std::vector<T>. The class is a `singleton', in that only
// one can exist, and access is provided via the static method
// G4SolidStore::GetInstance().
//
// All solids should be registered with G4SolidStore, and removed on their
// destruction. Intended principally for UI browser. The underlying
// container initially has a capacity of 100.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for std::vector<T>
//
// Member data:
//
// static G4SolidStore*
//   - Ptr to the single G4SolidStore

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 10.07.95 P.Kent  Initial version

#ifndef G4VSOLIDSTORE_HH
#define G4VSOLIDSTORE_HH

#include "g4std/vector"

#include "G4VSolid.hh"

class G4SolidStore : public G4std::vector<G4VSolid*>
{
  public:  // with description

    static void Register(G4VSolid* pSolid);
      // Add the solid to the collection.
    static void DeRegister(G4VSolid* pSolid);
      // Remove the solid from the collection.
    static G4SolidStore* GetInstance();
      // Get a ptr to the unique G4SolidStore, creating it if necessary.

    virtual ~G4SolidStore();
      // Default destructor.

  protected:

    G4SolidStore();

  private:

    static G4SolidStore* fgInstance;
};

#endif
