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
// 
// Class Description
//
// This class is a container for pointers to physics process objects.

// --------------------------------------------------------------------
#ifndef G4ProcessVector_hh
#define G4ProcessVector_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

class G4VProcess;

class G4ProcessVector 
{
  public:

    //  Constructors
    G4ProcessVector();
    explicit G4ProcessVector(size_t);
    G4ProcessVector(const G4ProcessVector &);

    //  Destructor.
    virtual ~G4ProcessVector();

    // assignment operator
    G4ProcessVector & operator=(const G4ProcessVector &right);
 
    // equality operator
    G4bool operator==(const G4ProcessVector &right) const;

  public: // With Description

    // Returns the number of items
    size_t entries() const;
    size_t length() const;
    size_t size() const;
    
    // Returns the position of the element
    size_t index(G4VProcess* aProcess) const;
 
    // Returns "true" if the element exists
    G4bool contains(G4VProcess* aProcess) const;

    // Insert an element
    G4bool insert(G4VProcess* aProcess);

    // Insert an element at i-th position
    G4bool insertAt(G4int i, G4VProcess* aProcess);
    
    // Remove and returns the i-th element
    G4VProcess* removeAt(G4int i);

    // Remove and returns the last element
    G4VProcess* removeLast();
    
    // Clear the collection by removing all items
    void clear();
    
    // returns const reference to the i-th item
    G4VProcess* const& operator[](G4int i) const;
    G4VProcess* const& operator()(G4int i) const;

    // returns  reference to the i-th item
    G4VProcess* & operator[](G4int i);
    G4VProcess* & operator()(G4int i);

  protected:

    typedef std::vector<G4VProcess*> G4ProcVector;

    G4ProcVector * pProcVector;
};

#include "G4ProcessVector.icc"

#endif
