// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneHandlerList.cc,v 1.1 2001-02-23 15:43:22 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4SceneHandlerList.hh"

G4SceneHandlerList::remove(G4VSceneHandler* sceneHandler) {
  G4SceneHandlerListIterator iSceneHandler;
  for (iSceneHandler = begin(); iSceneHandler != end(); ++iSceneHandler) {
    if (*iSceneHandler == sceneHandler) break;
  }
  if (iSceneHandler != end()) erase(iSceneHandler);
}
