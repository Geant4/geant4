//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: gtest01.cc 66241 2012-12-13 18:34:42Z gunter $
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

