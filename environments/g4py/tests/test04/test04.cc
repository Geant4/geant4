// $Id: test04.cc,v 1.2 2006-04-25 07:54:28 kmura Exp $
// ====================================================================
//   test04.cc
//
//   test of "Call Policies"
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "Step.hh"
#include "Track.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test04) {
  class_<Step>( "Step", "step class")
    .def(init<>())
    .add_property("x", &Step::GetX, &Step::SetX)
    ;

  class_<Track>( "Track", "track class")
    .def(init<>())
    .def("GetStep1", &Track::GetStep1,
	 return_internal_reference<>())
    .def("GetStep2", &Track::GetStep2,
	 return_value_policy<reference_existing_object>())
    // this is invalid, just for test
    .def("GetStep3", &Track::GetStep3,
	 return_value_policy<manage_new_object>())
    ;
}

