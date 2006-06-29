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
// $Id: pyG4UserStackingAction.cc,v 1.7 2006-06-29 15:31:44 gunter Exp $
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
      return f(boost::ref(aTrack));
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
  class_<CB_G4UserStackingAction, CB_G4UserStackingAction*, boost::noncopyable>
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

