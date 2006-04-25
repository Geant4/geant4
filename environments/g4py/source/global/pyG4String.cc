// $Id: pyG4String.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4String.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4String.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
std::ostream& operator<<(std::ostream& ostr, const G4String& astr)
{
  return ostr << astr.c_str();
}

// ====================================================================
// module definition
// ====================================================================
void export_G4String()
{
  class_<G4String>("G4String", "string class")
    .def(init<const G4String&>())
    .def(init<const char*>())
    .def(self_ns::str(self))
    .def(self + self)
    .def(self += self)
    .def(self += other<const char*>())
    .def(self == self)
    .def(self == other<const char*>())
    .def(self != self)
    .def(self != other<const char*>())
    ;
  
  implicitly_convertible<G4String, const char*>();
  implicitly_convertible<const char* ,G4String>();

  implicitly_convertible<G4String, std::string>();
  implicitly_convertible<std::string ,G4String>();
}
