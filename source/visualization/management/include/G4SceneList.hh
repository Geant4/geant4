// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.hh,v 1.1 1999-01-07 16:15:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4SCENELIST_HH
#define G4SCENELIST_HH

#include <rw/tpordvec.h>
#include "G4VScene.hh"

class G4SceneList: public RWTPtrOrderedVector<G4VScene> {
};

#endif
