// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VViewList.hh,v 2.0 1998/07/02 16:48:51 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  May 1996

#ifndef G4VVIEWLIST_HH
#define G4VVIEWLIST_HH

#include <rw/tpordvec.h>
#include "G4VView.hh"

class G4VViewList: public RWTPtrOrderedVector<G4VView> {};

#endif
