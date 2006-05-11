// $Id: gtest01.cc,v 1.1 2006-05-11 03:00:08 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   gtest01.cc
//
//   python wrapper for user application
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "QMaterials.hh"
#include "QDetectorConstruction.hh"
#include "QPhysicsList.hh"
#include "QPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "QEventAction.hh"

using namespace boost::python;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(gtest01) {
  class_<QMaterials>("QMaterials", "my material")
    .def("Construct", &QMaterials::Construct)
    ;

  class_<QDetectorConstruction, QDetectorConstruction*,
    bases<G4VUserDetectorConstruction> >
    ("QDetectorConstruction", "my detector")
    ;

  class_<QPhysicsList, QPhysicsList*,
    bases<G4VUserPhysicsList> >
    ("QPhysicsList", "my physics list")
    ;

  class_<QPrimaryGeneratorAction, QPrimaryGeneratorAction*,
    bases<G4VUserPrimaryGeneratorAction> >
    ("QPrimaryGeneratorAction", "my primary generator action")
    .def("GetParticleGun", &QPrimaryGeneratorAction::GetParticleGun,
         return_internal_reference<>())
    ;

  class_<QEventAction, QEventAction*,
    bases<G4UserEventAction> >
    ("QEventAction", "my event action")
    ;
}

