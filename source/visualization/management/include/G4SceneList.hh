// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.hh,v 1.3 1999-01-11 00:48:14 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  9th August 1998

#ifndef G4SCENELIST_HH
#define G4SCENELIST_HH

#include "G4Scene.hh"
#include <rw/tpordvec.h>

class G4SceneList: public RWTPtrOrderedVector <G4Scene> {};

#endif
