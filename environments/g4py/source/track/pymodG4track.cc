// $Id: pymodG4track.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4track.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Track();
void export_G4TrackStatus();
void export_G4Step();
void export_G4StepPoint();
void export_G4StepStatus();

BOOST_PYTHON_MODULE(G4track)
{
  export_G4Track();
  export_G4TrackStatus();
  export_G4Step();
  export_G4StepPoint();
  export_G4StepStatus();
}

