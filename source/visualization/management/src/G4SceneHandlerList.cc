// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneHandlerList.cc,v 1.2 2001-03-07 14:37:47 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4SceneHandlerList.hh"

void G4SceneHandlerList::remove(G4VSceneHandler* sceneHandler) {
  G4SceneHandlerListIterator iSceneHandler;
  for (iSceneHandler = begin(); iSceneHandler != end(); ++iSceneHandler) {
    if (*iSceneHandler == sceneHandler) break;
  }
  if (iSceneHandler != end()) erase(iSceneHandler);
}
