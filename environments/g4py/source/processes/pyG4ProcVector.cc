// $Id: pyG4ProcVector.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ProcVector.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4indexing.hh"
#include "G4VProcess.hh"

using namespace boost::python;

namespace pyG4ProcVector {

typedef std::vector<G4VProcess*> G4ProcVector;

};

using namespace pyG4ProcVector;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProcVector()
{
  class_<G4ProcVector> ("G4ProcVector", "process vector")
    .def(vector_indexing_suite<G4ProcVector>())
    ;
}

