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
// $Id: pyG4TwoVector.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4TwoVector.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 830
#include "G4TwoVector.hh"

using namespace boost::python;
using namespace CLHEP;

typedef G4TwoVector XXX; // ...

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4TwoVector {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_isNear, isNear, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_isParallel, isParallel, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_isOrthogonal, isOrthogonal, 1, 2)

}

using namespace pyG4TwoVector;

#endif

// ====================================================================
// module definition
// ====================================================================
void export_G4TwoVector()
{
#if G4VERSION_NUMBER >= 830

  class_<G4TwoVector>("G4TwoVector", "general 2-vector")
    // constructors
    .def(init<G4double>())
    .def(init<G4double, G4double>())
    .def(init<const XXX&>())

    // property
    .add_property("x", &XXX::x, &XXX::setX)
    .add_property("y", &XXX::y, &XXX::setY)

    // methods
    .def("set",      &XXX::set)
    .def("phi",      &XXX::phi)
    .def("mag",      &XXX::mag)
    .def("mag2",     &XXX::mag2)
    .def("r",        &XXX::r)
    .def("setPhi",   &XXX::setPhi)
    .def("setMag",   &XXX::setMag)
    .def("setR",     &XXX::setR)
    .def("setPolar", &XXX::setPolar)
    .def("howNear",  &XXX::howNear)
    .def("isNear",   &XXX::isNear,            f_isNear())
    .def("howParallel",   &XXX::howParallel)
    .def("isParallel",    &XXX::isParallel,   f_isParallel())
    .def("howOrthogonal", &XXX::howOrthogonal)
    .def("isOrthogonal",  &XXX::isOrthogonal, f_isOrthogonal())
    .def("unit",       &XXX::unit)
    .def("orthogonal", &XXX::orthogonal)
    .def("dot",        &XXX::dot)
    .def("angle",      &XXX::angle)
    .def("rotate",     &XXX::rotate)

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
    .def(self >  self)
    .def(self <  self)
    .def(self >= self)
    .def(self <= self)
    ;
#endif

}


