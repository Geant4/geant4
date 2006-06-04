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
// $Id: pyG4VPhysicsConstructor.cc,v 1.2 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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
				       SetPhysicsName, 0, 1);

};

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

