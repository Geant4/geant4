// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ViewerList.hh,v 1.2.2.1 1999/12/07 20:53:54 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// John Allison  May 1996

#ifndef G4VIEWERLIST_HH
#define G4VIEWERLIST_HH

#include "g4rw/tpordvec.h"
#include "G4VViewer.hh"

class G4ViewerList: public G4RWTPtrOrderedVector<G4VViewer> {};

#endif
