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
// $Id: pyG4VUserPrimaryGeneratorAction.cc,v 1.7 2006-06-29 15:35:33 gunter Exp $
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

