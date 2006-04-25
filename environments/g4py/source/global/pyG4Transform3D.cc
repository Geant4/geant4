// $Id: pyG4Transform3D.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Transform3D.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

using namespace boost::python;

typedef G4Transform3D XXX; // ...

// ====================================================================
// module definition
// ====================================================================
void export_G4Transform3D()
{
  class_<G4Transform3D>("G4Transform3D", "geometrical 3D transformation")
    // constructors
    .def(init<const G4RotationMatrix&, const G4ThreeVector&>())
    .def(init<const XXX&>())

    // property
    .add_property("xx", &XXX::xx)
    .add_property("xy", &XXX::xy)
    .add_property("xz", &XXX::xz)
    .add_property("yx", &XXX::yx)
    .add_property("yy", &XXX::yy)
    .add_property("yz", &XXX::yz)
    .add_property("zx", &XXX::zx)
    .add_property("zy", &XXX::zy)
    .add_property("zz", &XXX::zz)
    .add_property("dx", &XXX::dx)
    .add_property("dy", &XXX::dy)
    .add_property("dz", &XXX::dz)
    .def_readonly("Identity", &XXX::Identity)

    // methods
    .def("inverse",        &XXX::inverse)
    .def("getRotation" ,   &XXX::getRotation)
    .def("getTranslation", &XXX::getTranslation)

    // operators
    .def(self == self)
    .def(self != self)
    .def(self *  self)
    ;
}

