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
// $Id: pyG4PrimaryVertex.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetPrimary, GetPrimary, 0, 1);

};

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
