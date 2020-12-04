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
// ====================================================================
//   pyPhysListFactory.cc
//
//                                         2020 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4PhysListFactory.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_PhysListFactory()
{
  class_<G4PhysListFactory, G4PhysListFactory*>
    ("G4PhysListFactory", "phys list factory")
    .def("GetReferencePhysList", &G4PhysListFactory::GetReferencePhysList,
         return_internal_reference<>())
    .def("ReferencePhysList", &G4PhysListFactory::ReferencePhysList,
         return_internal_reference<>())
    .def("IsReferencePhysList", &G4PhysListFactory::IsReferencePhysList)
    .def("AvailablePhysLists", &G4PhysListFactory::AvailablePhysLists,
	 return_value_policy<reference_existing_object>())
    .def("AvailablePhysListsEM", &G4PhysListFactory::AvailablePhysListsEM,
	 return_value_policy<reference_existing_object>())
    ;
}
