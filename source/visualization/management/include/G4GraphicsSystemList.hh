// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GraphicsSystemList.hh,v 1.2 1999-01-09 16:30:38 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  2nd April 1996

#ifndef G4GRAPHICSSYSTEMLIST_HH
#define G4GRAPHICSSYSTEMLIST_HH

#include <rw/tpordvec.h>
#include "G4VGraphicsSystem.hh"

class G4GraphicsSystemList: public RWTPtrOrderedVector<G4VGraphicsSystem> {
};

#endif
