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
// $Id: G4VSolidRefArray.hh,v 1.1.4.1 2001/06/28 19:11:28 gunter Exp $
// GEANT4 tag $Name:  $
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
