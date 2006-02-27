// $Id: pymodG4interface.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4pyInterface.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================

void export_G4UImanager();
void export_G4UIterminal();

BOOST_PYTHON_MODULE(G4interface)
{
  export_G4UImanager();
  export_G4UIterminal();
}

