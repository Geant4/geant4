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
// G4OrderedTable
//
// Class description:
//
// Utility class, defining an ordered collection of vectors of <G4double>.

// Author: M.Maire (LAPP), September 1996
// Revisions: H.Kurashige (Kobe Univ.), January-September 2001
// --------------------------------------------------------------------
#ifndef G4OrderedTable_hh
#define G4OrderedTable_hh 1

#include <vector>

#include "G4DataVector.hh"
#include "globals.hh"

class G4OrderedTable : public std::vector<G4DataVector*>
{
 public:
  G4OrderedTable() = default;
  // Default constructor

  explicit G4OrderedTable(std::size_t cap);
  // Constructor given a 'capacity' defining the initial
  // number of elements (NULL pointers are filled up)

  virtual ~G4OrderedTable() = default;
  // Empty Destructor

  void clearAndDestroy();
  // Removes all elements and deletes all non-NULL pointers

  G4bool Store(const G4String& filename, G4bool ascii = false);
  // Stores OrderedTable in a file (returns false in case of failure)

  G4bool Retrieve(const G4String& filename, G4bool ascii = false);
  // Retrieves OrderedTable from a file (returns false in case of failure)

  friend std::ostream& operator<<(std::ostream& out, G4OrderedTable& table);
};

#endif
