// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GraphicsSystemList.hh,v 1.4 1999-12-15 14:54:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  2nd April 1996

#ifndef G4GRAPHICSSYSTEMLIST_HH
#define G4GRAPHICSSYSTEMLIST_HH

#include "g4rw/tpordvec.h"
#include "G4VGraphicsSystem.hh"

class G4GraphicsSystemList: public G4RWTPtrOrderedVector<G4VGraphicsSystem> {
};

#endif
