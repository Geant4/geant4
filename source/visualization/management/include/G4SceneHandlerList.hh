// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneHandlerList.hh,v 1.1 1999-01-09 16:30:42 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4SCENEHANDLERLIST_HH
#define G4SCENEHANDLERLIST_HH

#include <rw/tpordvec.h>
#include "G4VSceneHandler.hh"

class G4SceneHandlerList: public RWTPtrOrderedVector<G4VSceneHandler> {
};

#endif
