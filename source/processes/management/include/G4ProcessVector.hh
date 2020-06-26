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
// G4ProcessVector
//
// Class description:
//
// A container for pointers to physics process objects.

// Authors: G.Cosmo, H.Kurashige - 1998
// --------------------------------------------------------------------
#ifndef G4ProcessVector_hh
#define G4ProcessVector_hh 1

#include <vector>

#include "globals.hh"
#include "G4ios.hh"

class G4VProcess;

class G4ProcessVector 
{
  public:

    G4ProcessVector();
    explicit G4ProcessVector(std::size_t);
    G4ProcessVector(const G4ProcessVector &);
      // Constructors

    virtual ~G4ProcessVector();
      // Destructor

    G4ProcessVector& operator=(const G4ProcessVector& right);
      // Assignment operator
 
    inline G4bool operator==(const G4ProcessVector& right) const;
      // Equality operator

    inline std::size_t entries() const;
    inline std::size_t length() const;
    inline std::size_t size() const;
      // Returns the number of items
    
    std::size_t index(G4VProcess* aProcess) const;
      // Returns the position of the element
 
    G4bool contains(G4VProcess* aProcess) const;
      // Returns "true" if the element exists

    inline G4bool insert(G4VProcess* aProcess);
      // Insert an element

    G4bool insertAt(G4int i, G4VProcess* aProcess);
      // Insert an element at i-th position
    
    G4VProcess* removeAt(G4int i);
      // Remove and returns the i-th element

    inline G4VProcess* removeLast();
      // Remove and returns the last element
    
    inline void clear();
      // Clear the collection by removing all items
    
    inline G4VProcess* const& operator[](G4int i) const;
    inline G4VProcess* const& operator()(G4int i) const;
      // Returns const reference to the i-th item

    inline G4VProcess* & operator[](G4int i);
    inline G4VProcess* & operator()(G4int i);
      // Returns reference to the i-th item

  protected:

    using G4ProcVector = std::vector<G4VProcess*>;

    G4ProcVector* pProcVector = nullptr;
};

#include "G4ProcessVector.icc"

#endif
