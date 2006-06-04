//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: pydemo_wp.cc,v 1.2 2006-06-04 21:37:25 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pydemo_wp.cc
//
//   python wrapper for user application
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "MyMaterials.hh"
#include "MyDetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "G4VSensitiveDetector.hh"

using namespace boost::python;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(demo_wp) {
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

