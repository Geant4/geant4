// $Id: pyQmaterials.cc,v 1.1 2006-02-27 09:48:56 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyQmaterials.cc
//
//   [Qmaterials]
//   a site-module of Geant4Py
//
//   A sample set of materials
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

void Construct();

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(Qmaterials) {
  def("Construct", Construct);
}

