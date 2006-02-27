// $Id: pydemo_wp.cc,v 1.1 2006-02-27 09:44:35 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   gtest01.cc
//
//   python wrapper for user application
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "MyApplication.hh"
#include "MyMaterials.hh"
#include "MyDetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "G4VSensitiveDetector.hh"

using namespace boost::python;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(demo_wp) {
  class_<MyApplication>("MyApplication", "my application")
    .def("Configure", &MyApplication::Configure)
    ;

  class_<MyMaterials>("MyMaterials", "my material")
    .def("Construct", &MyMaterials::Construct)
    ;

  class_<MyDetectorConstruction, MyDetectorConstruction*,
    bases<G4VUserDetectorConstruction> >
    ("MyDetectorConstruction", "my detector")
    .def("SetSDtoScoreVoxel", &MyDetectorConstruction::SetSDtoScoreVoxel)
    ;

  class_<MyPhysicsList, MyPhysicsList*,
    bases<G4VUserPhysicsList> >
    ("MyPhysicsList", "my physics list")
    ;

}

