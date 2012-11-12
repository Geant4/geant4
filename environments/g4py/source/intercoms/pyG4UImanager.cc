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
// $Id: pyG4UImanager.cc,v 1.1 2006-08-08 05:20:57 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UImanager.cc
//
//   G4UImanager class is pure singleton, so it cannnot be exposed
//   from BPL. Functionality of G4UImanager is exposed in global
//   name space via wrappers.
//                                         2006 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4UImanager {

// ApplyCommand
G4int(G4UImanager::*f1_ApplyCommand)(const char*) = &G4UImanager::ApplyCommand;
G4int(G4UImanager::*f2_ApplyCommand)(const G4String&) = 
  &G4UImanager::ApplyCommand;

// CreateHTML
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_CreateHTML, CreateHTML, 0, 1);


//////////////////////////////////////////////
G4int ApplyUICommand_1(const G4String& cmdstr)
//////////////////////////////////////////////
{
  G4UImanager* UImgr= G4UImanager::GetUIpointer();
  G4int returnVal= UImgr-> ApplyCommand(cmdstr);
  if( returnVal == fCommandSucceeded ) return returnVal;

  G4int paramIndex= returnVal % 100;
  G4int commandStatus= returnVal - paramIndex;
 
  switch(commandStatus) {
    case fCommandSucceeded:
      break;

    case fCommandNotFound:
      G4cout << "command <" << UImgr-> SolveAlias(cmdstr)
	     << "> not found" << G4endl;
      break;

    case fIllegalApplicationState:
      G4cout << "illegal application state -- command refused" 
	     << G4endl;
      break;

    case fParameterOutOfRange:
      break;

    case fParameterOutOfCandidates:
      G4cout << "Parameter is out of candidate list (index "
	     << paramIndex << ")"
	     << G4endl;
      break;

    case fParameterUnreadable:
      G4cout << "Parameter is wrong type and/or is not omittable (index " 
	     << paramIndex << ")" << G4endl;
      break;

    case fAliasNotFound:
      break;

    default:
      G4cout << "command refused (" << commandStatus << ")" << G4endl;
      break;
  }

  return returnVal;
}

/////////////////////////////////////////////////
G4int ApplyUICommand_2(const std::string& cmdstr)
/////////////////////////////////////////////////
{
  return ApplyUICommand_1(cmdstr);
}

};

using namespace pyG4UImanager;

// ====================================================================
// module definition
// ====================================================================
void export_G4UImanager()
{
 class_<G4UImanager, boost::noncopyable>
   ("G4UImanager", "UI manager class", no_init)
   .def("GetUIpointer",  &G4UImanager::GetUIpointer,
        return_value_policy<reference_existing_object>())
   .staticmethod("GetUIpointer")
   // ---
   .def("GetCurrentValues", &G4UImanager::GetCurrentValues)
   .def("ExecuteMacroFile", &G4UImanager::ExecuteMacroFile)
   .def("ApplyCommand",     f1_ApplyCommand)
   .def("ApplyCommand",     f2_ApplyCommand)
   .def("CreateHTML",       &G4UImanager::CreateHTML, f_CreateHTML())
   .def("SetMacroSearchPath", &G4UImanager::SetMacroSearchPath)
   .def("GetMacroSearchPath", &G4UImanager::GetMacroSearchPath,
        return_value_policy<return_by_value>())
   // ---
   .def("SetPauseAtBeginOfEvent", &G4UImanager::SetPauseAtBeginOfEvent)
   .def("GetPauseAtBeginOfEvent", &G4UImanager::GetPauseAtBeginOfEvent)
   .def("SetPauseAtEndOfEvent",   &G4UImanager::SetPauseAtEndOfEvent)
   .def("GetPauseAtEndOfEvent",   &G4UImanager::GetPauseAtEndOfEvent)
   .def("SetVerboseLevel",        &G4UImanager::SetVerboseLevel)
   .def("GetVerboseLevel",        &G4UImanager::GetVerboseLevel)
   // ---
   .def("GetTree",   &G4UImanager::GetTree,
        return_value_policy<reference_existing_object>())
   ;

  // ---
  def("ApplyUICommand",    ApplyUICommand_1);
  def("ApplyUICommand",    ApplyUICommand_2);

}
