// $Id: test11.cc,v 1.2 2006-04-25 07:54:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test11.cc
//
//   test for the case for lack of default constructor
//
//                                         2005 Q
// ====================================================================
class AClass {
private:
  AClass() { }

public:
  AClass(int i)  { }
  AClass(const AClass& a) { }
  ~AClass() { }
};


// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;


BOOST_PYTHON_MODULE(test11)
{
  class_<AClass>("AClass", "a class", no_init)
    .def(init<int>())
    ;
}

