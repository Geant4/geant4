// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ViewerList.hh,v 1.1 1999-01-09 16:30:53 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4VIEWERLIST_HH
#define G4VIEWERLIST_HH

#include <rw/tpordvec.h>
#include "G4VViewer.hh"

class G4ViewerList: public RWTPtrOrderedVector<G4VViewer> {};

#endif
