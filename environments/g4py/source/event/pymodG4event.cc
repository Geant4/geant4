// $Id: pymodG4event.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4event.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4EventManager();
void export_G4StackManager();
void export_G4Event();
void export_G4UserEventAction();
void export_G4UserStackingAction();
void export_G4ClassificationOfNewTrack();
void export_G4ParticleGun();

BOOST_PYTHON_MODULE(G4event)
{
  export_G4EventManager();
  export_G4StackManager();
  export_G4Event();
  export_G4UserEventAction();
  export_G4UserStackingAction();
  export_G4ClassificationOfNewTrack();
  export_G4ParticleGun();
}

