// $Id: pyglobals.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyglobals.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4indexing.hh"
#include "globals.hh"
#include <vector>

using namespace boost::python;

namespace pyglobals {

typedef std::vector<G4int>    G4intVector;
typedef std::vector<G4double> G4doubleVector;
typedef std::vector<G4String> G4StringVector;

};

using namespace pyglobals;

// ====================================================================
// module definition
// ====================================================================
void export_globals()
{
  class_<G4intVector> ("G4intVector", "int vector")
    .def(vector_indexing_suite<G4intVector>())
    ;

  class_<G4doubleVector> ("G4doubleVector", "double vector")
    .def(vector_indexing_suite<G4doubleVector>())
    ;

  class_<G4StringVector> ("G4StringVector", "string vector")
    .def(vector_indexing_suite<G4StringVector>())
    ;
}

