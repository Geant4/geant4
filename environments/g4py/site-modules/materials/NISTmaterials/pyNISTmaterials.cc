// $Id: pyNISTmaterials.cc,v 1.2 2006-04-25 10:28:36 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyNISTmaterials.cc
//
//   [NISTmaterials]
//   a site-module of Geant4Py
//
//   NIST materials
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

void Construct();

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(NISTmaterials) {
  def("Construct", Construct);
}

