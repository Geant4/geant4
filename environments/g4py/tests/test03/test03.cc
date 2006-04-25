// $Id: test03.cc,v 1.2 2006-04-25 07:54:28 kmura Exp $
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

