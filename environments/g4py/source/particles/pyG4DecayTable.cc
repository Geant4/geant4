// $Id: pyG4DecayTable.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4DecayTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4DecayTable.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4DecayTable()
{
  class_<G4DecayTable, G4DecayTable*, boost::noncopyable>
    ("G4DecayTable", "decay table")
     // ---
     .def("DumpInfo",   &G4DecayTable::DumpInfo)
     ;

     // reduced functionality...
     // ...

}

