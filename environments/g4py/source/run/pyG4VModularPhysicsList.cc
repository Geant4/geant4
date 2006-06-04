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
// $Id: pyG4VModularPhysicsList.cc,v 1.2 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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
    .def("RegisterPhysis",     &G4VModularPhysicsList::RegisterPhysics)
    .def("GetPhysics",         f1_GetPhysics,
         return_value_policy<reference_existing_object>())
    .def("GetPhysics",         f2_GetPhysics,
         return_value_policy<reference_existing_object>())
    ;
}

