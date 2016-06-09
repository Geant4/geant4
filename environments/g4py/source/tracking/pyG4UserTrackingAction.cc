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
      f(boost::ref(atrack));
    } else {
      G4UserTrackingAction::PreUserTrackingAction(atrack);
    }
  }

  // PostUserTrackingAction
  void PostUserTrackingAction(const G4Track* atrack) {
    if(const override& f= get_override("PostUserTrackingAction")) {
      f(boost::ref(atrack));
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
  class_<CB_G4UserTrackingAction, CB_G4UserTrackingAction*, boost::noncopyable>
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

