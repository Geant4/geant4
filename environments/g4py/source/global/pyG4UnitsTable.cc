// $Id: pyG4UnitsTable.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UnitsTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UnitsTable.hh"

using namespace boost::python;


// ====================================================================
// module definition
// ====================================================================
void export_G4UnitsTable()
{
  class_<G4BestUnit>("G4BestUnit", "present best unit", no_init)
    .def(init<G4double, const G4String&>())
    .def(init<const G4ThreeVector&, const G4String&>())
    // ---
    //.def("GetValue",         & ...)
    .def("GetCategory",        &G4BestUnit::GetCategory,
	 return_internal_reference<>())
    .def("GetIndexOfCategory", &G4BestUnit::GetIndexOfCategory)
    .def(self_ns::str(self))
    ;

  implicitly_convertible<G4BestUnit, G4String>();

}

