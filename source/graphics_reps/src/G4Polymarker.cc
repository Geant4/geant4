// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polymarker.cc,v 1.1 1999-01-07 16:09:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  November 1996

#include "G4Polymarker.hh"

// Forward declaration for g++ linker.
#ifdef GNU_GCC
  #include <rw/tvvector.h>
  #include "G4Point3D.hh"
  template class RWTValVector<G4Point3D>;
#endif

G4Polymarker::G4Polymarker ():
fMarkerType (line)
{}

ostream& operator << (ostream& os, const G4Polymarker& marker) {
  os << "G4Polymarker: type: ";
  switch (marker.fMarkerType) {
  case G4Polymarker::line:
    os << "line"; break;
  case G4Polymarker::dots:
    os << "dots"; break;
  case G4Polymarker::circles:
    os << "circles"; break;
  case G4Polymarker::squares:
    os << "squares"; break;
  default:
    os << "unrecognised"; break;
  }
  os << "\n  ";
  os << (G4VMarker) marker;
  os << (G4Point3DList) marker;
  return os;
}
