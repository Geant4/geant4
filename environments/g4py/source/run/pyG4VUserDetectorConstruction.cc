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
// $Id: pyG4VUserDetectorConstruction.cc,v 1.6 2007-07-10 07:09:51 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VUserDetectorConstruction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VUserDetectorConstruction {

struct CB_G4VUserDetectorConstruction :
  G4VUserDetectorConstruction, wrapper<G4VUserDetectorConstruction> {

  G4VPhysicalVolume* Construct() {
    return get_override("Construct")();
  }
};

};

using namespace pyG4VUserDetectorConstruction;


// ====================================================================
// module definition
// ====================================================================
void export_G4VUserDetectorConstruction()
{
  class_<CB_G4VUserDetectorConstruction, CB_G4VUserDetectorConstruction*,
    boost::noncopyable>
    ("G4VUserDetectorConstruction",
     "base class of user detector construction")

    .def("Construct",
	 pure_virtual(&G4VUserDetectorConstruction::Construct),
         return_value_policy<reference_existing_object>())
    ;
}
