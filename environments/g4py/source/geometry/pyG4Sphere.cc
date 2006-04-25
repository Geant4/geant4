// $Id: pyG4Sphere.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Sphere.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Sphere.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Sphere()
{
  class_<G4Sphere, G4Sphere*, bases<G4VSolid> >
    ("G4Sphere", "Sphere solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
         G4double, G4double, G4double>())
    // ---
    .def("GetInsideRadius",     &G4Sphere::GetInsideRadius)
    .def("GetOuterRadius",      &G4Sphere::GetOuterRadius)
    .def("GetStartPhiAngle",    &G4Sphere::GetStartPhiAngle)
    .def("GetDeltaPhiAngle",    &G4Sphere::GetDeltaPhiAngle)
    .def("GetStartThetaAngle",  &G4Sphere::GetStartThetaAngle)
    .def("GetDeltaThetaAngle",  &G4Sphere::GetDeltaThetaAngle)
    .def("SetInsideRadius",     &G4Sphere::SetInsideRadius)
    .def("SetOuterRadius",      &G4Sphere::SetOuterRadius)
    .def("SetStartPhiAngle",    &G4Sphere::SetStartPhiAngle)
    .def("SetDeltaPhiAngle",    &G4Sphere::SetDeltaPhiAngle)
    .def("SetStartThetaAngle",  &G4Sphere::SetStartThetaAngle)
    .def("SetDeltaThetaAngle",  &G4Sphere::SetDeltaThetaAngle)
    .def("GetCubicVolume",      &G4Sphere::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
