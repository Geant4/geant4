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
// G4DataVector
//
// Class description:
//
// Utility class providing similar behaviour of std::vector<G4double>.
// It includes additional methods for compatibility with Rogue Wave
// collection.

// Author: H.Kurashige, 18 September 2001
// --------------------------------------------------------------------
#ifndef G4DataVector_hh
#define G4DataVector_hh 1

#include <fstream>
#include <iostream>
#include <vector>

#include "G4ios.hh"
#include "globals.hh"

class G4DataVector : public std::vector<G4double>
{
 public:
  G4DataVector() = default;
  // Default constructor.

  G4DataVector(const G4DataVector&) = default;
  G4DataVector(G4DataVector&&)      = default;
  // Default copy&move constructors.

  explicit G4DataVector(std::size_t cap);
  // Constructor given a 'capacity' defining the initial number of elements.

  G4DataVector(std::size_t cap, G4double value);
  // Constructor given a 'capacity' defining the initial number of elements
  // and initialising them to 'value'.

  virtual ~G4DataVector() = default;
  // Empty destructor

  G4DataVector& operator=(const G4DataVector&) = default;
  G4DataVector& operator=(G4DataVector&&) = default;
  // Default copy&move assignment operators.

  inline void insertAt(std::size_t, const G4double&);
  // Insert an element at given position

  inline std::size_t index(const G4double&) const;
  // Returns back index of the element same as given value

  inline G4bool contains(const G4double&) const;
  // Returns 'true' if it contains the element same as given value

  inline G4bool remove(const G4double&);
  // Removes the first element same as given value

  inline std::size_t removeAll(const G4double&);
  // Remove all elements same as given value

  enum
  {
    T_G4DataVector = 100
  };

  G4bool Store(std::ofstream& fOut, G4bool ascii = false);
  G4bool Retrieve(std::ifstream& fIn, G4bool ascii = false);
  // To store/retrieve persistent data to/from file streams.

  friend std::ostream& operator<<(std::ostream&, const G4DataVector&);
};

#include "G4DataVector.icc"

#endif
