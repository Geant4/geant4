// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneHandlerList.hh,v 1.3 1999-12-15 14:54:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4SCENEHANDLERLIST_HH
#define G4SCENEHANDLERLIST_HH

#include "g4rw/tpordvec.h"
#include "G4VSceneHandler.hh"

class G4SceneHandlerList: public G4RWTPtrOrderedVector<G4VSceneHandler> {
};

#endif
