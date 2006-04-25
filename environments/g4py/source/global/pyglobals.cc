// $Id: pyglobals.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyglobals.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4strstreambuf.hh"
#include "G4UImanager.hh"
#include "G4PyCoutDestination.hh"
#include "pyG4indexing.hh"
#include <vector>

using namespace boost::python;

extern G4strstreambuf G4coutbuf;
extern G4strstreambuf G4cerrbuf;

namespace pyglobals {

// G4cout/cerr are set to Python stdout
void SetG4PyCoutDestination()
{
  G4UImanager* UImgr= G4UImanager::GetUIpointer();
  G4PyCoutDestination* pycout= new G4PyCoutDestination();
  G4coutbuf.SetDestination(pycout);
  G4cerrbuf.SetDestination(pycout);
}

void ResetG4PyCoutDestination()
{
  G4UImanager* UImgr= G4UImanager::GetUIpointer();
  UImgr-> SetCoutDestination(0);
}

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
  def("SetG4PyCoutDestination",    SetG4PyCoutDestination);
  def("ResetG4PyCoutDestination",  ResetG4PyCoutDestination);

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

