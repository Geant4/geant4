// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DList.cc,v 1.1 1999-01-07 16:09:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

#include "G4Point3DList.hh"

#include "G4ios.hh"

ostream& operator << (ostream& os, const G4Point3DList& points)
{
  os << "G4Point3DList[" << points.entries() << "]: ";
  for (int i = 0; i < points.entries(); i++) os << points(i);
  return os;
}
