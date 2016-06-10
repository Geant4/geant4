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
// $Id: pyG4UserEventAction.cc 66892 2013-01-17 10:57:59Z gunter $
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
      f(boost::ref(anEvent));
    } else
      G4UserEventAction::BeginOfEventAction(anEvent);
  }

  // EndOfEventAction
  void EndOfEventAction(const G4Event* anEvent) {
    if(const override& f= get_override("EndOfEventAction")) {
      f(boost::ref(anEvent));
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
  class_<CB_G4UserEventAction, CB_G4UserEventAction*, boost::noncopyable>
    ( "G4UserEventAction", "event action class")
    
    .def("BeginOfEventAction", &G4UserEventAction::BeginOfEventAction,
	 &CB_G4UserEventAction::BeginOfEventAction)
    .def("EndOfEventAction", &G4UserEventAction::EndOfEventAction,
	 &CB_G4UserEventAction::EndOfEventAction)
    ;
}

