// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DList.cc,v 1.2 1999-05-12 16:11:03 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

#include "G4Point3DList.hh"

#include "G4ios.hh"

G4Point3DList::~G4Point3DList () {}

ostream& operator << (ostream& os, const G4Point3DList& points)
{
  os << "G4Point3DList[" << points.entries() << "]: ";
  for (int i = 0; i < points.entries(); i++) os << points(i);
  return os;
}
