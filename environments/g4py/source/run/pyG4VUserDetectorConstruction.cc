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
// $Id: pyG4VUserDetectorConstruction.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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
    get_override("Construct")();
  }
};

};

using namespace pyG4VUserDetectorConstruction;


// ====================================================================
// module definition
// ====================================================================
void export_G4VUserDetectorConstruction()
{
  class_<CB_G4VUserDetectorConstruction, boost::noncopyable>
    ("G4VUserDetectorConstruction",
     "base class of user detector construction")

    .def("Construct",
	 pure_virtual(&G4VUserDetectorConstruction::Construct),
         return_value_policy<reference_existing_object>())
    ;
}
