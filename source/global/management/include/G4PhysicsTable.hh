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
// $Id: G4PhysicsTable.hh 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
// Class description:
//
// G4PhysicsTable is an utility class for storage of pointers
// to G4PhysicsVector containers. It derives all functionalities
// of STL vector containers with the addition of few methods for
// compatibility with previous implementation based on Rogue-Wave
// pointer collections.
// The constructor given the 'capacity' of the table, pre-allocates
// memory for the specified value by invoking the STL's reserve()
// function, in order to avoid reallocation during insertions.
// G4PhysicsTable has a vector of boolean which are used 
// as 'recalc-needed' flags when processes calculate physics tables. 
// ------------------------------------------------------------
//
// History:
// -------
// - First implementation, based on object model of
//   2nd December 1995. G.Cosmo
// - 1st March 1996, modified. K.Amako
// - 24th February 2001, migration to STL vectors. H.Kurashige
// - 9th March 2001, added Store/RetrievePhysicsTable. H.Kurashige
// - 20th August 2004, added FlagArray and related methods   H.Kurashige
//-------------------------------------

#ifndef G4PhysicsTable_h
#define G4PhysicsTable_h 1

#include <vector>
#include "globals.hh"
#include "G4ios.hh"

class G4PhysicsVector;

class G4PhysicsTable : public std::vector<G4PhysicsVector*> 
{

  typedef std::vector<G4PhysicsVector*> G4PhysCollection;
  typedef std::vector<G4bool> G4FlagCollection;
 
 public: // with description

  G4PhysicsTable();
    // Default constructor.

  explicit G4PhysicsTable(size_t cap);
    // Constructor with capacity. Reserves memory for the
    // specified capacity.

  virtual ~G4PhysicsTable();
    // Destructor.
    // Does not invoke deletion of contained pointed collections.

  G4PhysicsVector*& operator()(size_t);
  G4PhysicsVector* const& operator()(size_t) const;
    // Access operators.

  void clearAndDestroy();
    // Removes all items and deletes them at the same time.

  void   push_back( G4PhysicsVector* );
  void   insert (G4PhysicsVector*);
    // Pushes new element to collection.

  void   insertAt (size_t, G4PhysicsVector*); 
    // insert element at the specified position in the collection.
  
  void   resize(size_t, G4PhysicsVector* vec = (G4PhysicsVector*)(0));
  // resize collection
 
  size_t entries() const;
  size_t length() const;
    // Return collection's size.

  G4bool isEmpty() const;
    // Flags if collection is empty or not.

  G4bool ExistPhysicsTable(const G4String& fileName) const;
    // Check if the specified file exists or not

  G4bool StorePhysicsTable(const G4String& filename, G4bool ascii=false);
    // Stores PhysicsTable in a file (returns false in case of failure).
  
  G4bool RetrievePhysicsTable(const G4String& filename, G4bool ascii=false);
    // Retrieves Physics from a file (returns false in case of failure).

  void ResetFlagArray();
    // Reset the array of flags and all flags are set "true" 
    // This flag is supposed to be used as "recalc-needed" flag
    //   associated with each physics vector  

  G4bool GetFlag(size_t i) const;
  void   ClearFlag(size_t i);
    // Get/Clear the flag for the 'i-th' physics vector    
   
  friend std::ostream& operator<<(std::ostream& out, G4PhysicsTable& table);

 protected:

  G4PhysicsVector* CreatePhysicsVector(G4int type);  
  G4FlagCollection vecFlag; 

 private:

  G4PhysicsTable(const G4PhysicsTable&) = delete;
  G4PhysicsTable& operator=(const G4PhysicsTable&) = delete;
    // Private copy constructor and assignment operator.

};

typedef G4PhysicsTable::iterator G4PhysicsTableIterator;

#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.icc"

#endif
