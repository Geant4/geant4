#include <iostream.h>
#include "FallibleBase.hh"

void
FallibleBase::ErrUsedInInvalidState::writeMessage(ostream& os) const {
  os << "in function Fallible::operator T()" << endl;
  os << "Object used in invalid state";
}


