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
// $Id: pyG4VModularPhysicsList.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4VModularPhysicsList.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VModularPhysicsList.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VModularPhysicsList {

struct CB_G4VModularPhysicsList :
  G4VModularPhysicsList, wrapper<G4VModularPhysicsList> {

  void SetCuts() {
    get_override("SetCuts")();
  }

};

// GetPhysics()
const G4VPhysicsConstructor*
      (G4VModularPhysicsList::*f1_GetPhysics)(G4int) const
  = &G4VModularPhysicsList::GetPhysics;
const G4VPhysicsConstructor*
      (G4VModularPhysicsList::*f2_GetPhysics)(const G4String&) const
  = &G4VModularPhysicsList::GetPhysics;

}

using namespace pyG4VModularPhysicsList;

// ====================================================================
// module definition
// ====================================================================
void export_G4VModularPhysicsList()
{
  class_<CB_G4VModularPhysicsList, bases<G4VUserPhysicsList>,
    boost::noncopyable>
    ("G4VModularPhysicsList", "base class of modular physics list")
    // ---
    .def("SetCuts",            pure_virtual(&G4VModularPhysicsList::SetCuts))
    .def("ConstructParticle",  &G4VModularPhysicsList::ConstructParticle)
    .def("ConstructProcess",   &G4VModularPhysicsList::ConstructProcess)
    // ---
    .def("RegisterPhysics",     &G4VModularPhysicsList::RegisterPhysics)
    .def("GetPhysics",         f1_GetPhysics,
         return_value_policy<reference_existing_object>())
    .def("GetPhysics",         f2_GetPhysics,
         return_value_policy<reference_existing_object>())
    ;
}

