// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DList.cc,v 1.4 2001-05-03 11:06:10 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

#include "G4Point3DList.hh"

#include "G4ios.hh"

G4Point3DList::~G4Point3DList () {}

G4std::ostream& operator << (G4std::ostream& os, const G4Point3DList& points)
{
  os << "G4Point3DList[" << points.size() << "]: ";
  for (size_t i = 0; i < points.size(); i++) os << points[i];
  return os;
}
