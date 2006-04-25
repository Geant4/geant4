// $Id: pyG4RotationMatrix.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4RotationMatrix.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RotationMatrix.hh"

using namespace boost::python;

typedef G4RotationMatrix XXX; // ...

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4RotationMatrix {

XXX&(XXX::*f1_rotate)(G4double, const G4ThreeVector&)= &XXX::rotate;
XXX&(XXX::*f2_rotate)(G4double, const G4ThreeVector*)= &XXX::rotate;

};

using namespace pyG4RotationMatrix;

// ====================================================================
// module definition
// ====================================================================
void export_G4RotationMatrix()
{
  class_<G4RotationMatrix>("G4RotationMatrix", "rotation matrix")
    // constructors
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
    .def_readonly("IDENTITY", &XXX::IDENTITY)

    // methods
    .def("colX",       &XXX::colX)
    .def("colY",       &XXX::colY)
    .def("colZ",       &XXX::colZ)
    .def("rowX",       &XXX::rowX)
    .def("rowY",       &XXX::rowY)
    .def("rowZ",       &XXX::rowZ)
    .def("getPhi",     &XXX::getPhi)
    .def("getTheta",   &XXX::getTheta)
    .def("getPsi",     &XXX::getPsi)
    .def("phi",        &XXX::phi)
    .def("theta",      &XXX::theta)
    .def("psi",        &XXX::psi)
    .def("getDelta",   &XXX::getDelta)
    .def("getAxis",    &XXX::getAxis)
    .def("delta",      &XXX::axis)
    .def("axis",       &XXX::delta)
    .def("phiX",       &XXX::phiX)
    .def("phiY",       &XXX::phiY)
    .def("phiZ",       &XXX::phiZ)
    .def("thetaX",     &XXX::thetaX)
    .def("thetaY",     &XXX::thetaY)
    .def("thetaZ",     &XXX::thetaZ)
    .def("setPhi",     &XXX::setPhi)
    .def("setTheta",   &XXX::setTheta)
    .def("setPsi",     &XXX::setPsi)
    .def("setAxis",    &XXX::setAxis)
    .def("setDelta",   &XXX::setDelta)
    .def("isIdentity", &XXX::isIdentity)
    .def("rotateX",    &XXX::rotateX,
	 return_value_policy<reference_existing_object>())
    .def("rotateY",    &XXX::rotateY,
	 return_value_policy<reference_existing_object>())
    .def("rotateZ",    &XXX::rotateZ,
	 return_value_policy<reference_existing_object>())
    .def("rotate",     f1_rotate,
	 return_value_policy<reference_existing_object>())
    .def("rotate",     f2_rotate,
	 return_value_policy<reference_existing_object>())
    .def("rotateAxes", &XXX::rotateAxes,
	 return_value_policy<reference_existing_object>())
    .def("inverse",     &XXX::inverse)
    .def("invert",     &XXX::invert,
	 return_value_policy<reference_existing_object>())

    // operators
    .def(self_ns::str(self))
    .def(self == self)
    .def(self != self)
    .def(self >  self)
    .def(self <  self)
    .def(self >= self)
    .def(self <= self)
    .def(self *  self)
    .def(self *  G4ThreeVector())
    .def(self *= self)
    ;
}

