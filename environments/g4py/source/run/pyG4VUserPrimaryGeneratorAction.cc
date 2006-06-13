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
// $Id: pyG4VUserPrimaryGeneratorAction.cc,v 1.5 2006-06-13 09:54:43 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VUserPrimaryGeneratorAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VUserPrimaryGeneratorAction {

struct CB_G4VUserPrimaryGeneratorAction :
  G4VUserPrimaryGeneratorAction, wrapper<G4VUserPrimaryGeneratorAction> {
  
  void GeneratePrimaries(G4Event* anEvent) {
    //get_override("GeneratePrimaries")(&anEvent);
    get_override("GeneratePrimaries")(boost::ref(anEvent));
  }
};

};

using namespace pyG4VUserPrimaryGeneratorAction;

// ====================================================================
// module definition
// ====================================================================
void export_G4VUserPrimaryGeneratorAction()
{
  class_<CB_G4VUserPrimaryGeneratorAction, CB_G4VUserPrimaryGeneratorAction*,
    boost::noncopyable> 
    ("G4VUserPrimaryGeneratorAction", 
     "base class of user primary generator action")

    .def("GeneratePrimaries", 
	 pure_virtual(&G4VUserPrimaryGeneratorAction::GeneratePrimaries))
    ;
}

