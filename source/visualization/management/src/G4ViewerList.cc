// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ViewerList.cc,v 1.2 2001-03-07 14:37:50 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4ViewerList.hh"

void G4ViewerList::remove(G4VViewer* viewer) {
  G4ViewerListIterator iViewer;
  for (iViewer = begin(); iViewer != end(); ++iViewer) {
    if (*iViewer == viewer) break;
  }
  if (iViewer != end()) erase(iViewer);
}
