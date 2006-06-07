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
// $Id: pyG4UserEventAction.cc,v 1.4 2006-06-07 05:22:05 kmura Exp $
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
  class_<CB_G4UserEventAction, CB_G4UserEventAction*, boost::noncopyable>
    ( "G4UserEventAction", "event action class")
    
    .def("BeginOfEventAction", &G4UserEventAction::BeginOfEventAction,
	 &CB_G4UserEventAction::BeginOfEventAction)
    .def("EndOfEventAction", &G4UserEventAction::EndOfEventAction,
	 &CB_G4UserEventAction::EndOfEventAction)
    ;
}

