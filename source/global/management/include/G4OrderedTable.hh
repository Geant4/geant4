// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OrderedTable.hh,v 1.1 1999-01-07 16:09:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//
//	This class is setting up an ordered collection of 
//      ordered vectors of <G4double>
//	30 September 1996, M.Maire
//
// ------------------------------------------------------------

#ifndef G4OrderedTable_h
#define G4OrderedTable_h 1

#include "globals.hh"
#include <rw/tvordvec.h>
#include <rw/tpordvec.h>

class G4ValVector : public RWTValOrderedVector<G4double>
{

  public:

      G4ValVector(size_t capac=RWDEFAULT_CAPACITY)
        : RWTValOrderedVector<G4double>(capac) {;}

      virtual ~G4ValVector() {;}


      G4bool operator==(const G4ValVector &right) const
      {
        return (this == (G4ValVector *) &right);
      }

};

typedef RWTPtrOrderedVector<G4ValVector> G4OrderedTable;

#endif
