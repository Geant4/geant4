// $Id: pyG4ElementTable.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ElementTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4indexing.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ElementTable()
{
  class_<G4ElementTable> ("G4ElementTable", "element table")
    .def(vector_indexing_suite<G4ElementTable>())
    .def(self_ns::str(self))
    ;
}

