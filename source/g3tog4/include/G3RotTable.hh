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
// $Id: G3RotTable.hh,v 1.12 2001-07-11 09:58:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// G3 rotations table.
// Maps G3 rotations indices to their G4RotationMatrix object counterparts.

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3ROTTABLEH_HH
#define G3ROTTABLEH_HH 1

#include "G3RotTableEntry.hh"

#include "globals.hh"

#include "g4rw/tpordvec.h"

class G4Material;

typedef G4RWTPtrOrderedVector<G3RotTableEntry>  G3RotMatrixVector;

class G3RotTable
{
  public:  // with description

    G3RotTable();
    virtual ~G3RotTable();
    
    // methods

    G4RotationMatrix* Get(G4int id) const;
    void Put(G4int id, G4RotationMatrix* matrix);
    void Clear();

  private:

    G3RotMatrixVector*  fRotVector;
};

extern G3RotTable G3Rot;

#endif
