// $Id: pyG4UserLimits.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UserLimits.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserLimits.hh"
#include "G4Track.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4UserLimits()
{
  class_<G4UserLimits, G4UserLimits*>
    ("G4UserLimits", "user step limitations")
    // ---
    .def(init<G4double>())
    .def(init<G4double, G4double>())
    .def(init<G4double, G4double, G4double>())
    .def(init<G4double, G4double, G4double, G4double>())
    .def(init<G4double, G4double, G4double, G4double, G4double>())
    // ---
    .def(init<const G4String&>())
    .def(init<const G4String&, G4double>())
    .def(init<const G4String&, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double>())
    // ---
    .def("GetUserMaxTrackLength", &G4UserLimits::GetUserMaxTrackLength)
    .def("GetUserMaxTime",        &G4UserLimits::GetUserMaxTime)
    .def("GetUserMinEkine",       &G4UserLimits::GetUserMinEkine)
    .def("GetUserMinRange",       &G4UserLimits::GetUserMinRange)
    // ---
    .def("SetMaxAllowedStep",     &G4UserLimits::SetMaxAllowedStep)
    .def("SetUserMaxTrackLength", &G4UserLimits::SetUserMaxTrackLength)
    .def("SetUserMaxTime",        &G4UserLimits::SetUserMaxTime)
    .def("SetUserMinEkine",       &G4UserLimits::SetUserMinEkine)
    .def("SetUserMinRange",       &G4UserLimits::SetUserMinRange)
    // ---
    .def("GetType",               &G4UserLimits::GetType,
         return_internal_reference<>())
    .def("SetType",               &G4UserLimits::SetType)
    ;
}
