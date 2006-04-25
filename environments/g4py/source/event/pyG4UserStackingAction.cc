// $Id: pyG4UserStackingAction.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UserStackingAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserStackingAction.hh"
#include "G4Track.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4UserStackingAction {

struct CB_G4UserStackingAction : G4UserStackingAction,
				 wrapper<G4UserStackingAction> {
  
  // ClassifyNewTrack
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack) {
    if(const override& f= get_override("ClassifyNewTrack")) {
      return f(aTrack);
    } else
      return G4UserStackingAction::ClassifyNewTrack(aTrack);
  }

  // NewState
  void NewStage() {
    if(const override& f= get_override("NewStage")) {
      f();
    } else
      G4UserStackingAction::NewStage();
  }

  // PrepareNewEvent
  void PrepareNewEvent() {
    if(const override& f= get_override("PrepareNewEvent")) {
      f();
    } else
      G4UserStackingAction::PrepareNewEvent();
  }

};

};

using namespace pyG4UserStackingAction;

// ====================================================================
// module definition
// ====================================================================
void export_G4UserStackingAction()
{
  class_<G4UserStackingAction, boost::noncopyable>
    ("G4UserStackingAction", "stacking action class")
    // ---
    .def("ClassifyNewTrack",  &G4UserStackingAction::ClassifyNewTrack,
	 &CB_G4UserStackingAction::ClassifyNewTrack)
    .def("NewStage",          &G4UserStackingAction::NewStage,
	 &CB_G4UserStackingAction::NewStage)
    .def("PrepareNewEvent",   &G4UserStackingAction::PrepareNewEvent,
	 &CB_G4UserStackingAction::PrepareNewEvent)
    ;
}

