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
// $Id: G4OrderedTable.hh,v 1.9 2001-09-17 08:17:57 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
// Sep. 1996  : M.Maire 
// Jan. 2001  : H.Kurashige
//              - G4ValVector is replaced with G4DataVector 
//              - Migrated to G4std::vector<G4DataVector*>.
// Sep. 2001  : H.Kurashige
//              - Add
//                 G4bool Store(const G4String&, G4bool)
//                 G4bool Retrieve(const G4String&, G4bool);
//                 ostream& operator<<(ostream&, G4OrderedTable&)
//
// Class Description:
//
//	Utility class, defining an ordered collection of vectors
//      of <G4double>.

// ------------------------------------------------------------

#ifndef G4OrderedTable_h
#define G4OrderedTable_h 1

#include "globals.hh"
#include "g4std/vector"
class G4DataVector;

class G4OrderedTable : public G4std::vector<G4DataVector*> 
{

 public: // with description

  G4OrderedTable();
    // Deafult constructor.

  G4OrderedTable(size_t capacity);
    // Constructor given a 'capacity' defining the initial
    // number of elements (NULL pointers are filled up)

  virtual ~G4OrderedTable(){;}
    // Empty Destructor

  inline void clearAndDestroy();
    // Removes all elements and deletes all non-NULL pointers

  G4bool Store(const G4String& filename, G4bool ascii=false);
    // Stores OrderedTable in a file (returns false in case of failure).
  
  G4bool Retrieve(const G4String& filename, G4bool ascii=false);
    // Retrieves OrderedTable from a file (returns false in case of failure).

  friend G4std::ostream& operator<<(G4std::ostream& out, G4OrderedTable& table);

};

typedef G4OrderedTable::iterator G4OrderedTableIterator;

#include "G4DataVector.hh"

inline
G4OrderedTable::G4OrderedTable()
  : G4std::vector<G4DataVector*>()
{
}

inline
G4OrderedTable::G4OrderedTable(size_t capacity)
  : G4std::vector<G4DataVector*>(capacity, (G4DataVector*)(0) )
{
}

inline 
void G4OrderedTable::clearAndDestroy() 
{
  G4DataVector* a;
  while (size()>0) {
    a = back();
    pop_back();
    for (iterator i=begin(); i!=end(); i++){
      if (*i==a) {
	erase(i);
	i--;
      }
    } 
    if ( a )  delete a;    
  } 
}


#endif


