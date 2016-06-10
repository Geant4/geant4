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
// $Id: pyG4VPhysicsConstructor.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4VPhysicsConstructor.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VPhysicsConstructor.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VPhysicsConstructor {

class CB_G4VPhysicsConstructor :
    public G4VPhysicsConstructor,
    public wrapper<G4VPhysicsConstructor> {

public:
  CB_G4VPhysicsConstructor(): G4VPhysicsConstructor() { }
  CB_G4VPhysicsConstructor(const G4String& name)
    : G4VPhysicsConstructor(name) { }

  void ConstructParticle() {
    get_override("ConstructParticle")();
  }

  void ConstructProcess() {
    get_override("ConstructProcess")();
  }

};

// SetPhysicsName()
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetPhysicsName,
				                               SetPhysicsName, 0, 1)

}

using namespace pyG4VPhysicsConstructor;

// ====================================================================
// module definition
// ====================================================================
void export_G4VPhysicsConstructor()
{
  class_<CB_G4VPhysicsConstructor, boost::noncopyable>
    ("G4VPhysicsConstructor",
     "base class of user physics constructor")
    // ----
    .def(init<const G4String&>())
    // ---
    .def("ConstructParticle",
         pure_virtual(&G4VPhysicsConstructor::ConstructParticle))
    .def("ConstructProcess",
         pure_virtual(&G4VPhysicsConstructor::ConstructProcess))
    // ---
    .def("SetPhysicsName", &G4VPhysicsConstructor::SetPhysicsName,
	 f_SetPhysicsName())
    .def("GetPhysicsName", &G4VPhysicsConstructor::GetPhysicsName,
         return_value_policy<return_by_value>())
    .def("SetVerboseLevel", &G4VPhysicsConstructor::SetVerboseLevel)
    .def("GetVerboseLevel", &G4VPhysicsConstructor::GetVerboseLevel)
    ;
}

