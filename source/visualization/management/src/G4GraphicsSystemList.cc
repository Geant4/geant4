// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GraphicsSystemList.cc,v 1.2 2001-03-07 14:37:45 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4GraphicsSystemList.hh"

void G4GraphicsSystemList::remove(G4VGraphicsSystem* graphicsSystem) {
  G4GraphicsSystemListIterator iGS;
  for (iGS = begin(); iGS != end(); ++iGS) {
    if (*iGS == graphicsSystem) break;
  }
  if (iGS != end()) erase(iGS);
}
