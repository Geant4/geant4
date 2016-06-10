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
// $Id: pyG4UIcommand.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4UIcommand.cc
//
//                                         2006 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UIcommand.hh"
#include "G4UImessenger.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4UIcommand {

// GetStateList (returning python list)
list f_GetStateList(G4UIcommand* acommand)
{
  list pyStateList;
  std::vector<G4ApplicationState>* stateList= acommand->GetStateList();

  for( size_t i=0; i< stateList->size(); i++) {
    pyStateList.append(&(*stateList)[i]);
  }

  return pyStateList;
}

}

using namespace pyG4UIcommand;

// ====================================================================
// module definition
// ====================================================================
void export_G4UIcommand()
{
  class_<G4UIcommand, G4UIcommand*>
    ("G4UIcommand", "UI command")
    .def(init<const char*, G4UImessenger*>())
    // ---
    .def("GetCurrentValue",     &G4UIcommand::GetCurrentValue)
    .def("IsAvailable",         &G4UIcommand::IsAvailable)
    .def("List",                &G4UIcommand::List)
    .def("GetRange",            &G4UIcommand::GetRange,
         return_value_policy<return_by_value>())
    .def("GetGuidanceEntries",  &G4UIcommand::GetGuidanceEntries)
    .def("GetGuidanceLine",     &G4UIcommand::GetGuidanceLine,
         return_value_policy<return_by_value>())
    .def("GetCommandPath",      &G4UIcommand::GetCommandPath,
         return_value_policy<return_by_value>())
    .def("GetCommandName",      &G4UIcommand::GetCommandName,
         return_value_policy<return_by_value>())
    .def("GetParameterEntries", &G4UIcommand::GetParameterEntries)
    .def("GetParameter",        &G4UIcommand::GetParameter,
         return_value_policy<reference_existing_object>())
    .def("GetStateList",        f_GetStateList)
    .def("GetTitle",            &G4UIcommand::GetTitle)
    ;
}
