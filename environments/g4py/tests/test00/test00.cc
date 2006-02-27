// $Id: test00.cc,v 1.1 2006-02-27 10:05:24 kmura Exp $
// ====================================================================
//   test00.cc
//   Hallo World of Boost.Python
//                                         2005 Q
// ====================================================================

// C function to be wrappred
char const* greet()
{
   return "Hallo World!!";
}


// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test00)
{
  def("greet", greet);
}

