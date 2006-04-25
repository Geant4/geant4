// $Id: pyG4UserEventAction.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UserEventAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserEventAction.hh"
#include "G4Event.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
struct CB_G4UserEventAction : G4UserEventAction, 
			     wrapper<G4UserEventAction> {
  // BeginOfEventAction
  void BeginOfEventAction(const G4Event* anEvent) {
    if(const override& f= get_override("BeginOfEventAction")) {
      f(anEvent);
    } else
      G4UserEventAction::BeginOfEventAction(anEvent);
  }

  // EndOfEventAction
  void EndOfEventAction(const G4Event* anEvent) {
    if(const override& f= get_override("EndOfEventAction")) {
      f(anEvent);
    } else {
      G4UserEventAction::EndOfEventAction(anEvent);
    }    
  }
};


// ====================================================================
// module definition
// ====================================================================
void export_G4UserEventAction()
{
  class_<CB_G4UserEventAction, boost::noncopyable>
    ( "G4UserEventAction", "event action class")
    
    .def("BeginOfEventAction", &G4UserEventAction::BeginOfEventAction,
	 &CB_G4UserEventAction::BeginOfEventAction)
    .def("EndOfEventAction", &G4UserEventAction::EndOfEventAction,
	 &CB_G4UserEventAction::EndOfEventAction)
    ;
}

