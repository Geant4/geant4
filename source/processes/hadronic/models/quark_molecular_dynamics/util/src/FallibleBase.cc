#include "g4std/iostream"
#include "FallibleBase.hh"

void
FallibleBase::ErrUsedInInvalidState::writeMessage(G4std::ostream& os) const {
  os << "in function Fallible::operator T()" << G4endl;
  os << "Object used in invalid state";
}


