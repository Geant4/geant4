// $Id: pymodG4digits_hits.cc,v 1.2 2006-04-25 08:09:44 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4digits_hits.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================

void export_G4VSensitiveDetector();

BOOST_PYTHON_MODULE(G4digits_hits)
{
  export_G4VSensitiveDetector();
}

