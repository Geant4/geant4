// $Id: pymodG4tracking.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4tracking.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4TrackingManager();
void export_G4UserSteppingAction();
void export_G4UserTrackingAction();

BOOST_PYTHON_MODULE(G4tracking)
{
  export_G4TrackingManager();
  export_G4UserSteppingAction();
  export_G4UserTrackingAction();
}

