// $Id: pymodG4materials.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4materials.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4Version.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Material();
void export_G4MaterialTable();
void export_G4Element();
void export_G4ElementTable();
void export_G4NistManager();

BOOST_PYTHON_MODULE(G4materials)
{
  export_G4Material();
  export_G4MaterialTable();
  export_G4Element();
  export_G4ElementTable();
#if G4VERSION_NUMBER >= 710
  export_G4NistManager();
#endif
}

