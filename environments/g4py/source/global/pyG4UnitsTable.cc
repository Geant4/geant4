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
// $Id: pyG4UnitsTable.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UnitsTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UnitsTable.hh"
#include "pyG4indexing.hh"

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

