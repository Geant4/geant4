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
// $Id: pyG4Cons.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Cons.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Cons.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Cons()
{
  class_<G4Cons, G4Cons*, bases<G4VSolid> >
    ("G4Cons", "Cone solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double, G4double, G4double>())
    // ---
    .def("GetInnerRadiusMinusZ", &G4Cons::GetInnerRadiusMinusZ)
    .def("GetOuterRadiusMinusZ", &G4Cons::GetOuterRadiusMinusZ)
    .def("GetInnerRadiusPlusZ",  &G4Cons::GetInnerRadiusPlusZ)
    .def("GetOuterRadiusPlusZ",  &G4Cons::GetOuterRadiusPlusZ)
    .def("GetZHalfLength",       &G4Cons::GetZHalfLength)
    .def("GetStartPhiAngle",     &G4Cons::GetStartPhiAngle)
    .def("GetDeltaPhiAngle",     &G4Cons::GetDeltaPhiAngle)
    .def("SetInnerRadiusMinusZ", &G4Cons::SetInnerRadiusMinusZ)
    .def("SetOuterRadiusMinusZ", &G4Cons::SetOuterRadiusMinusZ)
    .def("SetInnerRadiusPlusZ",  &G4Cons::SetInnerRadiusPlusZ)
    .def("SetOuterRadiusPlusZ",  &G4Cons::SetOuterRadiusPlusZ)
    .def("SetZHalfLength",       &G4Cons::SetZHalfLength)
    .def("SetStartPhiAngle",     &G4Cons::SetStartPhiAngle)
    .def("SetDeltaPhiAngle",     &G4Cons::SetDeltaPhiAngle)
    .def("GetCubicVolume",       &G4Cons::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

