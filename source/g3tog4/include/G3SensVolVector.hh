// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3SensVolVector.hh,v 1.2 1999-12-05 17:50:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// vector of logical volumes that were defined with
// tracking medium with ISVOL=1
//
// by I.Hrivnacova, 27 Sep 99

#include "g4rw/tpordvec.h"
#include "G4LogicalVolume.hh"

typedef G4RWTPtrOrderedVector<G4LogicalVolume> G3SensVolVector;

extern G3SensVolVector G3SensVol;
