// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogVolRefArray.hh,v 1.1 2000/11/17 05:10:02 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

// class description:
//
//	This is a class which implements G4RWTPtrVector with
//      G4LogicalVolume.
//

// History:
// 00.11.17 Y.Morita  Initial version

#ifndef G4LogVolRefArray_hh
#define G4LogVolRefArray_hh 1

#include "g4rw/tpvector.h"
#include "globals.hh"

#include "G4LogicalVolume.hh"

typedef G4RWTPtrVector<G4LogicalVolume> G4LogVolRefVArray;

class G4LogVolRefArray 
{
  public: // with description
      G4LogVolRefArray();
      ~G4LogVolRefArray();
      //  The constructor and the destructor.

  private:
      G4LogVolRefVArray  transLogVolPtrs;

  public: // with description
      inline void Insert(G4int i, G4LogicalVolume* aVol)
      { transLogVolPtrs[i] = aVol; }
      // Insert aVol as i-th element.
      inline G4LogicalVolume* Get(G4int i)
      { return transLogVolPtrs[i]; }
      // return the i-th element.
      inline void Resize(size_t n)
      { transLogVolPtrs.resize(n); }
      // Resize the vector with n.

};

G4LogVolRefArray::G4LogVolRefArray()
{;}

G4LogVolRefArray::~G4LogVolRefArray()
{;}

#endif
