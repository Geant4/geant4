// $Id: pyG4VUserDetectorConstruction.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
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
