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
// $Id: pyG4ThreeVector.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4ThreeVector.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

using namespace boost::python;
using namespace CLHEP;

typedef G4ThreeVector XXX; // ...

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ThreeVector {

G4double(XXX::*f1_theta)() const= &XXX::theta;
G4double(XXX::*f2_theta)(const XXX&) const = &XXX::theta;

G4double(XXX::*f1_cosTheta)() const= &XXX::cosTheta;
G4double(XXX::*f2_cosTheta)(const XXX&) const = &XXX::cosTheta;

G4double(XXX::*f1_cos2Theta)() const= &XXX::cos2Theta;
G4double(XXX::*f2_cos2Theta)(const XXX&) const = &XXX::cos2Theta;

G4double(XXX::*f1_perp2)() const= &XXX::perp2;
G4double(XXX::*f2_perp2)(const XXX&) const = &XXX::perp2;

G4double(XXX::*f1_perp)() const= &XXX::perp;
G4double(XXX::*f2_perp)(const XXX&) const = &XXX::perp;

G4double(XXX::*f1_angle)() const= &XXX::angle;
G4double(XXX::*f2_angle)(const XXX&) const = &XXX::angle;

G4double(XXX::*f1_eta)() const= &XXX::eta;
G4double(XXX::*f2_eta)(const XXX&) const = &XXX::eta;

XXX(XXX::*f1_project)() const= &XXX::project;
XXX(XXX::*f2_project)(const XXX&) const = &XXX::project;

XXX(XXX::*f1_perpPart)() const= &XXX::perpPart;
XXX(XXX::*f2_perpPart)(const XXX&) const = &XXX::perpPart;

G4double(XXX::*f1_rapidity)() const= &XXX::rapidity;
G4double(XXX::*f2_rapidity)(const XXX&) const = &XXX::rapidity;

G4double(XXX::*f1_polarAngle)(const XXX&) const= &XXX::polarAngle;
G4double(XXX::*f2_polarAngle)(const XXX&, const XXX&) const = &XXX::polarAngle;

G4double(XXX::*f1_azimAngle)(const XXX&) const= &XXX::azimAngle;
G4double(XXX::*f2_azimAngle)(const XXX&, const XXX&) const = &XXX::azimAngle;

XXX&(XXX::*f1_rotate)(G4double, const XXX&)= &XXX::rotate;
XXX&(XXX::*f2_rotate)(const XXX&, G4double)= &XXX::rotate;
XXX&(XXX::*f3_rotate)(const HepAxisAngle&)= &XXX::rotate;
XXX&(XXX::*f4_rotate)(const HepEulerAngles&)= &XXX::rotate;
XXX&(XXX::*f5_rotate)(G4double, G4double, G4double)= &XXX::rotate;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_isNear, isNear, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_isParallel, isParallel, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_isOrthogonal, isOrthogonal, 1, 2)

}

using namespace pyG4ThreeVector;

// ====================================================================
// module definition
// ====================================================================
void export_G4ThreeVector()
{
  class_<G4ThreeVector>("G4ThreeVector", "general 3-vector")
    // constructors
    .def(init<G4double>())
    .def(init<G4double, G4double>())
    .def(init<G4double, G4double, G4double>())
    .def(init<const XXX&>())

    // property
    .add_property("x", &XXX::x, &XXX::setX)
    .add_property("y", &XXX::y, &XXX::setY)
    .add_property("z", &XXX::z, &XXX::setZ)

    // methods
    .def("set",      &XXX::set)
    .def("phi",      &XXX::phi)
    .def("mag",      &XXX::mag)
    .def("mag2",     &XXX::mag2)
    .def("setPhi",   &XXX::setPhi)
    .def("setTheta", &XXX::setTheta)
    .def("setMag",   &XXX::setMag)
    .def("setPerp",  &XXX::setPerp)
    .def("setCylTheta", &XXX::setCylTheta)
    .def("howNear",  &XXX::howNear)
    .def("deltaR",   &XXX::deltaR)
    .def("unit",     &XXX::unit)
    .def("orthogonal", &XXX::orthogonal)
    .def("dot",      &XXX::dot)
    .def("cross",    &XXX::cross)
    .def("pseudoRapidity", &XXX::pseudoRapidity)
    .def("setEta",   &XXX::setEta)
    .def("setCylEta",&XXX::setCylEta)
    .def("setRThetaPhi", &XXX::setRThetaPhi)
    .def("setREtaPhi",   &XXX::setREtaPhi)
    .def("setRhoPhiZ",   &XXX::setRhoPhiZ)
    .def("setRhoPhiEta", &XXX::setRhoPhiEta)
    .def("getX",     &XXX::getX)
    .def("getY",     &XXX::getY)
    .def("getZ",     &XXX::getZ)
    .def("getR",     &XXX::getR)
    .def("getTheta", &XXX::getTheta)
    .def("getPhi",   &XXX::getPhi)
    .def("r",        &XXX::r)
    .def("rho",      &XXX::rho)
    .def("getRho",   &XXX::getRho)
    .def("getEta",   &XXX::getEta)
    .def("setR",     &XXX::setR)
    .def("setRho",   &XXX::setRho)
    .def("compare",  &XXX::compare)
    .def("diff2",    &XXX::diff2)
    .def("setTolerance",  &XXX::setTolerance)
    .staticmethod("setTolerance")
    .def("getTolerance",  &XXX::getTolerance)
    .staticmethod("getTolerance")
    .def("isNear",       &XXX::isNear,       f_isNear())
    .def("isParallel",   &XXX::isParallel,   f_isParallel())
    .def("isOrthogonal", &XXX::isOrthogonal, f_isOrthogonal())
    .def("howParallel",   &XXX::howParallel)
    .def("howOrthogonal", &XXX::howOrthogonal)
    .def("beta",     &XXX::beta)
    .def("gamma",    &XXX::gamma)
    .def("deltaPhi", &XXX::deltaPhi)
    .def("coLinearRapidity", &XXX::coLinearRapidity)
    .def("theta",     f1_theta)
    .def("theta",     f2_theta)
    .def("cosTheta",  f1_cosTheta)
    .def("cosTheta",  f2_cosTheta)
    .def("cos2Theta", f1_cos2Theta)
    .def("cos2Theta", f2_cos2Theta)
    .def("perp2",     f1_perp2)
    .def("perp2",     f2_perp2)
    .def("angle",     f1_angle)
    .def("angle",     f2_angle)
    .def("eta",       f1_eta)
    .def("eta",       f2_eta)
    .def("project",   f1_project)
    .def("project",   f2_project)
    .def("perpPart",  f1_perpPart)
    .def("perpPart",  f2_perpPart)
    .def("rapidity",  f1_rapidity)
    .def("rapidity",  f2_rapidity)
    .def("polarAngle",f1_polarAngle)
    .def("polarAngle",f2_polarAngle)
    .def("azimAngle", f1_azimAngle)
    .def("azimAngle", f2_azimAngle)
    .def("rotateX",   &XXX::rotateX,
	 return_value_policy<reference_existing_object>())
    .def("rotateY",   &XXX::rotateY,
	 return_value_policy<reference_existing_object>())
    .def("rotateZ",   &XXX::rotateZ,
	 return_value_policy<reference_existing_object>())
    .def("rotateUz", &XXX::rotateUz,
	 return_value_policy<reference_existing_object>())
    .def("transform",&XXX::transform,
	 return_value_policy<reference_existing_object>())
    .def("rotate",   f1_rotate,
	 return_value_policy<reference_existing_object>())
    .def("rotate",   f2_rotate,
	 return_value_policy<reference_existing_object>())
    .def("rotate",   f5_rotate,
	 return_value_policy<reference_existing_object>())

    // operators
    .def(self_ns::str(self))
    .def(self == self)
    .def(self != self)
    .def(self += self)
    .def(self -= self)
    .def(self -  self)
    .def(self + self)
    .def(self * self)
    .def(self * G4double())
    .def(self / G4double())
    .def(G4double() * self)
    .def(self *= G4double())
    .def(self /= G4double())
    .def(self >  self)
    .def(self <  self)
    .def(self >= self)
    .def(self <= self)
    ;
}

