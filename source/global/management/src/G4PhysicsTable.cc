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
  if (this != &right) {
    G4std::vector<G4PhysicsVector*>::const_iterator itr;
    for (itr=right.begin(); itr!=right.end(); ++itr){
      G4std::vector<G4PhysicsVector*>::push_back(*itr);
    }
  }
  return *this;
}

G4PhysicsTable::~G4PhysicsTable()
{
  clear();
}
