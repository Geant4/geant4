// $Id: test03.cc,v 1.1 2006-02-27 10:05:24 kmura Exp $
// ====================================================================
//   test03.cc
//
//   Singleton class (w/o public constructor) can be exposed?
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "AClass.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test03) {
  class_<AClass>("AClass", "Singleton")
    .def("GetPointer", &AClass::GetPointer,
	 return_value_policy<manage_new_object>())
    ;
}

