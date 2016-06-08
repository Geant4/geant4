// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSolidRefArray.hh,v 1.1 2000/11/17 05:10:03 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

// class description:
//
//	This is a class which implements G4RWTPtrVector with
//      G4VSolid.
//

// History:
// 00.11.17 Y.Morita  Initial version

#ifndef G4VSolidRefArray_hh
#define G4VSolidRefArray_hh 1

#include "g4rw/tpvector.h"
#include "globals.hh"

#include "G4VSolid.hh"

typedef G4RWTPtrVector<G4VSolid> G4VSolidRefVArray;

class G4VSolidRefArray 
{
  public: // with description
      G4VSolidRefArray();
      ~G4VSolidRefArray();
      //  The constructor and the destructor.

  private:
      G4VSolidRefVArray  transSolidPtrs;

  public: // with description
      inline void Insert(G4int i, G4VSolid* aSolid)
      { transSolidPtrs[i] = aSolid; }
      // Insert aVol as i-th element.
      inline G4VSolid* Get(G4int i)
      { return transSolidPtrs[i]; }
      // return the i-th element.
      inline void Resize(size_t n)
      { transSolidPtrs.resize(n); }
      // Resize the vector with n.

};

G4VSolidRefArray::G4VSolidRefArray()
{;}

G4VSolidRefArray::~G4VSolidRefArray()
{;}

#endif
