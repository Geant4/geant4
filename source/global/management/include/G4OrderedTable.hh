// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OrderedTable.hh,v 1.6 2001-02-05 10:37:12 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
// Sep. 1996  : M.Maire 
// Jan. 2001  : H.Kurashige
//              - G4ValVector is replaced with G4DataVector 
//              - Migrated to G4std::vector<G4DataVector*> from
//                G4RWTPtrOrderedVector
//
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
#include "G4DataVector.hh"

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

};

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
