// $Id: pymodG4run.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4run.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================

void export_G4RunManager();
void export_G4RunManagerKernel();
void export_G4Run();
void export_G4UserRunAction();
void export_G4VUserPrimaryGeneratorAction();
void export_G4VUserDetectorConstruction();
void export_G4VUserPhysicsList();
void export_G4VModularPhysicsList();
void export_G4VPhysicsConstructor();

BOOST_PYTHON_MODULE(G4run) 
{
  export_G4RunManager();
  export_G4RunManagerKernel();
  export_G4Run();
  export_G4UserRunAction();
  export_G4VUserPrimaryGeneratorAction();
  export_G4VUserDetectorConstruction();
  export_G4VUserPhysicsList();
  export_G4VModularPhysicsList();
  export_G4VPhysicsConstructor();
}
