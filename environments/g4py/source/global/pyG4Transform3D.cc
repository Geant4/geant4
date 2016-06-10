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
// $Id: pyG4Transform3D.cc 66892 2013-01-17 10:57:59Z gunter $
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

