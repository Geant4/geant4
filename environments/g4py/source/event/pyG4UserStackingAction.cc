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
// $Id: pyG4UserStackingAction.cc,v 1.4 2006-06-07 05:22:05 kmura Exp $
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
  class_<G4UserStackingAction, G4UserStackingAction*, boost::noncopyable>
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

