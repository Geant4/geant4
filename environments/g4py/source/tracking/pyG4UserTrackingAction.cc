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
// $ID$
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UserTrackingAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserTrackingAction.hh"
#include "G4Track.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
struct CB_G4UserTrackingAction : G4UserTrackingAction,
                                 wrapper<G4UserTrackingAction> {

  // PreUserTrackingAction
  void PreUserTrackingAction(const G4Track* atrack) {
    if(const override& f= get_override("PreUserTrackingAction")) {
      f(atrack);
    } else {
      G4UserTrackingAction::PreUserTrackingAction(atrack);
    }
  }

  // PostUserTrackingAction
  void PostUserTrackingAction(const G4Track* atrack) {
    if(const override& f= get_override("PostUserTrackingAction")) {
      f(atrack);
    } else {
      G4UserTrackingAction::PostUserTrackingAction(atrack);
    }
  }

};


// ====================================================================
// module definition
// ====================================================================
void export_G4UserTrackingAction()
{
  class_<CB_G4UserTrackingAction, boost::noncopyable>
    ("G4UserTrackingAction", "tracking action class")
    // ---
    .def("PreUserTrackingAction", 
	 &G4UserTrackingAction::PreUserTrackingAction,
         &CB_G4UserTrackingAction::PreUserTrackingAction)
    .def("PostUserTrackingAction", 
	 &G4UserTrackingAction::PostUserTrackingAction,
         &CB_G4UserTrackingAction::PostUserTrackingAction)
    ;
}

