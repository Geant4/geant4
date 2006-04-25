// $Id: pyG4UserRunAction.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
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

