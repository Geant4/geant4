// $Id: pyG4UniformMagField.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UniformMagField.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UniformMagField.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4UniformMagField()
{
  class_<G4UniformMagField, G4UniformMagField*, 
    bases<G4Field, G4MagneticField> >
    ("G4UniformMagField", "uniform magnetic field", no_init)
    // constructors
    .def(init<const G4ThreeVector&>())
    .def(init<const G4double, G4double, G4double>())
    // ---
    .def("SetFieldValue",         &G4UniformMagField::SetFieldValue)
    .def("GetConstantFieldValue", &G4UniformMagField::GetConstantFieldValue)
    ;
}
