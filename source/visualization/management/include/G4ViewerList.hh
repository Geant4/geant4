// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ViewerList.hh,v 1.4 2001-02-23 15:43:19 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4VIEWERLIST_HH
#define G4VIEWERLIST_HH

#include "g4std/vector"
#include "G4VViewer.hh"

class G4ViewerList: public G4std::vector<G4VViewer*> {
public:
  remove(G4VViewer*);
};

typedef G4std::vector<G4VViewer*>::iterator G4ViewerListIterator;
typedef G4std::vector<G4VViewer*>::const_iterator G4ViewerListConstIterator;

#endif
