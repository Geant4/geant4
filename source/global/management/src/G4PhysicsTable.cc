// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsTable.cc,v 1.2 2001-03-05 09:31:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation
//
//	G4PhysicsTable
//
// ------------------------------------------------------------

#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.hh"

G4PhysicsTable::G4PhysicsTable()
  : G4std::vector<G4PhysicsVector*>()
{
}

G4PhysicsTable::G4PhysicsTable(size_t capacity)
  : G4std::vector<G4PhysicsVector*>()
{
  reserve(capacity);
}

G4PhysicsTable::G4PhysicsTable(const G4PhysicsTable& right)
{
  *this = right;
}

G4PhysicsTable& G4PhysicsTable::operator=(const G4PhysicsTable& right)
{
  if (this != &right)
  {
    G4std::vector<G4PhysicsVector*>::const_iterator itr;
    for (itr=right.begin(); itr!=right.end(); ++itr)
    {
      G4std::vector<G4PhysicsVector*>::push_back(*itr);
    }
  }
  return *this;
}

G4PhysicsTable::~G4PhysicsTable()
{
  clear();
}
