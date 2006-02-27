// $Id: test08.cc,v 1.1 2006-02-27 10:05:25 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test08.cc
//
//   test for static member function
//
//                                         2005 Q
// ====================================================================

class AClass {
public:
  AClass() { }
  ~AClass() { }
  static int AMethod() { return 1; } 

};

// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test08)
{
  class_<AClass>( "AClass", "a class")
    .def(init<>())
    .def("AMethod", &AClass::AMethod)
    .staticmethod("AMethod")
    ;
}

