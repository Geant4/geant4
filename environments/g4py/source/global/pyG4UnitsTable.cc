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
// $Id: pyG4UnitsTable.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4UnitsTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4indexing.hh"
#include "G4UnitsTable.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4UnitsTable()
{
  class_<G4UnitsTable>("G4UnitsTable", "Units Table")
    .def(vector_indexing_suite<G4UnitsTable>())
    ;

  class_<G4UnitsContainer>("G4UnitsContainer", "Units Container")
    .def(vector_indexing_suite<G4UnitsContainer>())
    ;

  class_<G4UnitDefinition, boost::noncopyable>
    ("G4UnitDefinition", "Unit Definition", no_init)
    .def(init<const G4String&, const G4String&, const G4String&, G4double>())
    // ---
    .def("GetName",           &G4UnitDefinition::GetName,
	 return_value_policy<return_by_value>())
    .def("GetSymbol",         &G4UnitDefinition::GetSymbol,
	 return_value_policy<return_by_value>())
    .def("GetValue",          &G4UnitDefinition::GetValue)
    .def("PrintDefinition",   &G4UnitDefinition::PrintDefinition)
    // ---
    .def("BuildUnitsTable",   &G4UnitDefinition::BuildUnitsTable)
    .staticmethod("BuildUnitsTable")
    .def("PrintUnitsTable",   &G4UnitDefinition::PrintUnitsTable)
    .staticmethod("PrintUnitsTable")
    .def("GetUnitsTable",   &G4UnitDefinition::GetUnitsTable,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetUnitsTable")
    // ---
    .def("GetValueOf",          &G4UnitDefinition::GetValueOf)
    .staticmethod("GetValueOf")
    .def("GetCategory",         &G4UnitDefinition::GetCategory)
    .staticmethod("GetCategory")    
    ;

  class_<G4UnitsCategory, boost::noncopyable>
    ("G4UnitsCategory", "Units Category", no_init)
     .def(init<const G4String&>())
     // ---
    .def("GetName",           &G4UnitsCategory::GetName,
	 return_value_policy<return_by_value>())
    .def("GetUnitsList",      &G4UnitsCategory::GetUnitsList,
         return_value_policy<reference_existing_object>())
    .def("GetNameMxLen",      &G4UnitsCategory::GetNameMxLen)
    .def("GetSymbMxLen",      &G4UnitsCategory::GetSymbMxLen)
    .def("UpdateNameMxLen",   &G4UnitsCategory::UpdateNameMxLen)
    .def("UpdateSymbMxLen",   &G4UnitsCategory::UpdateSymbMxLen)
    .def("PrintCategory",     &G4UnitsCategory::PrintCategory)
    ;

  class_<G4BestUnit>("G4BestUnit", "present best unit", no_init)
    .def(init<G4double, const G4String&>())
    .def(init<const G4ThreeVector&, const G4String&>())
    // ---
    .def("GetCategory",        &G4BestUnit::GetCategory,
    	 return_value_policy<return_by_value>())
    .def("GetIndexOfCategory", &G4BestUnit::GetIndexOfCategory)
    .def(self_ns::str(self))
    ;

}

