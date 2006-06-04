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
// $Id: gtest01.cc,v 1.2 2006-06-04 21:35:59 kmura Exp $
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

