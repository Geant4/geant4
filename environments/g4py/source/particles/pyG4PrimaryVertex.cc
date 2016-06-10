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
// $Id: pyG4PrimaryVertex.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4PrimaryVertex.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4PrimaryVertex.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4PrimaryVertex {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetPrimary, GetPrimary, 0, 1)

}

using namespace pyG4PrimaryVertex;

// ====================================================================
// module definition
// ====================================================================
void export_G4PrimaryVertex()
{
  class_<G4PrimaryVertex, G4PrimaryVertex*>
    ("G4PrimaryVertex", "primary vertex")
    // ---
    .add_property("X0", &G4PrimaryVertex::GetX0)
    .add_property("Y0", &G4PrimaryVertex::GetY0)
    .add_property("Z0", &G4PrimaryVertex::GetZ0)
    .add_property("T0", &G4PrimaryVertex::GetT0)
    // ---
    .def("GetPosition", &G4PrimaryVertex::GetPosition,
         return_value_policy<return_by_value>())
    .def("GetX0",       &G4PrimaryVertex::GetX0)
    .def("GetY0",       &G4PrimaryVertex::GetY0)
    .def("GetZ0",       &G4PrimaryVertex::GetZ0)
    .def("GetT0",       &G4PrimaryVertex::GetT0)
    .def("GetNumberOfParticle", &G4PrimaryVertex::GetNumberOfParticle)
    .def("GetPrimary",  &G4PrimaryVertex::GetPrimary,
	 return_internal_reference<>(), f_GetPrimary())
    .def("GetWeight",   &G4PrimaryVertex::GetWeight)
    .def("SetWeight",   &G4PrimaryVertex::SetWeight)
    .def("Print", &G4PrimaryVertex::Print)
     ;
}
