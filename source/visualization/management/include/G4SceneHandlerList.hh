// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneHandlerList.hh,v 1.4 2001-02-23 15:43:16 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4SCENEHANDLERLIST_HH
#define G4SCENEHANDLERLIST_HH

#include "g4std/vector"
#include "G4VSceneHandler.hh"

class G4SceneHandlerList: public G4std::vector<G4VSceneHandler*> {
public:
  remove(G4VSceneHandler*);
};

typedef G4std::vector<G4VSceneHandler*>::iterator G4SceneHandlerListIterator;
typedef G4std::vector<G4VSceneHandler*>::const_iterator
        G4SceneHandlerListConstIterator;

#endif
