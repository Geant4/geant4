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
// $Id: pyG4UserRunAction.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4UserRunAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserRunAction.hh"
#include "G4Run.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
struct CB_G4UserRunAction : G4UserRunAction, wrapper<G4UserRunAction> {
  // BeginOfRunAction
  void BeginOfRunAction(const G4Run* aRun) {
    if(const override& f= get_override("BeginOfRunAction")) {
      f(boost::ref(aRun));
    } else
      G4UserRunAction::BeginOfRunAction(aRun);
  }

  // EndOfRunAction
  void EndOfRunAction(const G4Run* aRun) {
    if(const override& f= get_override("EndOfRunAction")) {
      f(boost::ref(aRun));
    } else {
      G4UserRunAction::EndOfRunAction(aRun);
    }
  }
};


// ====================================================================
// module definition
// ====================================================================
void export_G4UserRunAction()
{
  class_<CB_G4UserRunAction, CB_G4UserRunAction*, boost::noncopyable>
    ( "G4UserRunAction", "run action class")
    // ---
    .def("BeginOfRunAction", &G4UserRunAction::BeginOfRunAction,
         &CB_G4UserRunAction::BeginOfRunAction)
    .def("EndOfRunAction", &G4UserRunAction::EndOfRunAction,
         &CB_G4UserRunAction::EndOfRunAction)

    // reduced functionality...
    //.def("GenerateRun",  &G4UserRunAction::GenerateRun) // virtual
    ;
}

