// $Id: pyG4Step.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Step.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Step.hh"

using namespace boost::python;


// ====================================================================
// module definition
// ====================================================================
void export_G4Step()
{
  class_<G4Step, G4Step*>("G4Step", "step class")
    // ---
    .def("GetTrack",                 &G4Step::GetTrack,
         return_value_policy<reference_existing_object>())
    .def("GetPreStepPoint",          &G4Step::GetPreStepPoint,
         return_internal_reference<>())
    .def("GetPostStepPoint",         &G4Step::GetPostStepPoint,
         return_internal_reference<>())
    .def("GetTotalEnergyDeposit",    &G4Step::GetTotalEnergyDeposit)
    .def("GetStepLength",            &G4Step::GetStepLength)
    .def("GetDeltaPosition",         &G4Step::GetDeltaPosition)
    .def("GetDeltaTime",             &G4Step::GetDeltaTime)
    .def("GetDeltaMomentum",         &G4Step::GetDeltaMomentum)
    .def("GetDeltaEnergy",           &G4Step::GetDeltaEnergy)
    ;
}

