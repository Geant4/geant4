// $Id: test07.cc,v 1.2 2006-04-25 07:54:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test07.cc
//
//   test for emuns
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

enum Color { red=100, blue, yellow };


// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test07)
{
  enum_<Color>("Color")
    .value("red", red)
    .value("blue", blue)
    .value("yellow", yellow)
    ;
}

