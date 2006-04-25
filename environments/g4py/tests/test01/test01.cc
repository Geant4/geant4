// $Id: test01.cc,v 1.2 2006-04-25 07:54:28 kmura Exp $
// ====================================================================
//   test01.cc
//
//   exposure of a simple class
//   - constructor
//   - multiple constructors
//   - default value in function arguments
//   - setter/getter
//   - method
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "AClass.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test01) {
  class_<AClass>( "AClass", "a simple class")
    .def(init<>())
    .def(init<int>())
    .def(init<int, double>())
    .add_property("ival", &AClass::GetIVal, &AClass::SetIVal)
    .add_property("dval", &AClass::GetDVal, &AClass::SetDVal)
    .def("AMethod", &AClass::AMethod)
    ;
}

