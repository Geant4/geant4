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
// $Id: pyG4Transform3D.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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

