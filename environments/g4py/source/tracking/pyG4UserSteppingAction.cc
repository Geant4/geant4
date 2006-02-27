// $ID$
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UserSteppingAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
struct CB_G4UserSteppingAction : G4UserSteppingAction,
				 wrapper<G4UserSteppingAction> {
  // UserSteppingAction
  void UserSteppingAction(const G4Step* astep) {
    if(const override& f= get_override("UserSteppingAction")) {
      f(astep);
    } else {
      G4UserSteppingAction::UserSteppingAction(astep);
    }
  }
};


// ====================================================================
// module definition
// ====================================================================
void export_G4UserSteppingAction()
{
  class_<CB_G4UserSteppingAction, boost::noncopyable>
    ("G4UserSteppingAction", "stepping action class")

    .def("UserSteppingAction", &G4UserSteppingAction::UserSteppingAction,
         &CB_G4UserSteppingAction::UserSteppingAction)
    ;
}

