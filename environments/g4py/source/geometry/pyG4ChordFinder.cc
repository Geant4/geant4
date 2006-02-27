// $Id: pyG4ChordFinder.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ChordFinder.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ChordFinder.hh"
#include "G4MagneticField.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ChordFinder {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetVerbose, SetVerbose, 0, 1);

};

using namespace pyG4ChordFinder;

// ====================================================================
// module definition
// ====================================================================
void export_G4ChordFinder()
{
  class_<G4ChordFinder, G4ChordFinder*, boost::noncopyable>
    ("G4ChordFinder", "chord finder class", no_init)
    // constructor
    .def(init<G4MagInt_Driver*>())
    .def(init<G4MagneticField*>())
    .def(init<G4MagneticField*, G4double>())
    .def(init<G4MagneticField*, G4double, G4MagIntegratorStepper*>())
    // ---
    .def("GetDeltaChord",   &G4ChordFinder::GetDeltaChord)
    .def("SetDeltaChord",   &G4ChordFinder::SetDeltaChord)
    // ---
    .def("PrintStatistics", &G4ChordFinder::PrintStatistics)
    .def("SetVerbose",      &G4ChordFinder::SetVerbose, f_SetVerbose())
    ;
}
