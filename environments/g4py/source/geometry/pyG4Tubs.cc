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
// $Id: pyG4Tubs.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Tubs.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Tubs.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Tubs()
{
  class_<G4Tubs, G4Tubs*, bases<G4VSolid> >
    ("G4Tubs", "Tube solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double>())
    // ---
    .def("GetInnerRadius",   &G4Tubs::GetInnerRadius)
    .def("GetOuterRadius",   &G4Tubs::GetOuterRadius)
    .def("GetZHalfLength",   &G4Tubs::GetZHalfLength)
    .def("GetStartPhiAngle", &G4Tubs::GetStartPhiAngle)
    .def("GetDeltaPhiAngle", &G4Tubs::GetDeltaPhiAngle)
    .def("SetInnerRadius",   &G4Tubs::SetInnerRadius)
    .def("SetOuterRadius",   &G4Tubs::SetOuterRadius)
    .def("SetZHalfLength",   &G4Tubs::SetZHalfLength)
    .def("SetStartPhiAngle", &G4Tubs::SetStartPhiAngle)
    .def("SetDeltaPhiAngle", &G4Tubs::SetDeltaPhiAngle)
    .def("GetRMin",          &G4Tubs::GetRMin)
    .def("GetRMax",          &G4Tubs::GetRMax)
    .def("GetDz",            &G4Tubs::GetDz)
    .def("GetSPhi",          &G4Tubs::GetSPhi)
    .def("GetDPhi",          &G4Tubs::GetDPhi)
    .def("GetCubicVolume",   &G4Tubs::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

