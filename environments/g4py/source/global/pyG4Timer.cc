// $Id: pyG4Timer.cc,v 1.1 2006-04-25 08:13:51 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Timer.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Timer.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Timer()
{
  class_<G4Timer>("G4Timer", "Timer")
    // ---
    .def("Start",            &G4Timer::Start)
    .def("Stop",             &G4Timer::Stop)
    .def("IsValid",          &G4Timer::IsValid)
    .def("GetRealElapsed",   &G4Timer::GetRealElapsed)
    .def("GetSystemElapsed", &G4Timer::GetSystemElapsed)
    .def("GetUserElapsed",   &G4Timer::GetUserElapsed)
    ;
}
