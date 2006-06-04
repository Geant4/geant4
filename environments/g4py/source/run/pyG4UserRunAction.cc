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
// $Id: pyG4UserRunAction.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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
      f(aRun);
    } else
      G4UserRunAction::BeginOfRunAction(aRun);
  }

  // EndOfRunAction
  void EndOfRunAction(const G4Run* aRun) {
    if(const override& f= get_override("EndOfRunAction")) {
      f(aRun);
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
  class_<CB_G4UserRunAction, boost::noncopyable>
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

