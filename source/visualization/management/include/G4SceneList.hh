// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.hh,v 1.5 1999-12-15 14:54:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  9th August 1998

#ifndef G4SCENELIST_HH
#define G4SCENELIST_HH

#include "G4Scene.hh"
#include "g4rw/tpordvec.h"

class G4SceneList: public G4RWTPtrOrderedVector <G4Scene> {};

#endif
