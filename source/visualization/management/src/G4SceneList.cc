// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.cc,v 1.4 2001-03-07 14:37:48 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4SceneList.hh"

void G4SceneList::remove(G4Scene* scene) {
  G4SceneListIterator iScene;
  for (iScene = begin(); iScene != end(); ++iScene) {
    if (*iScene == scene) break;
  }
  if (iScene != end()) erase(iScene);
}
