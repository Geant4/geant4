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
// $Id: pyG4Trap.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Trap.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Trap.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Trap()
{
  class_<G4Trap, G4Trap*, bases<G4VSolid> >
    ("G4Trap", "Generic trapezoild soild class", no_init)
    // constructors
    .def(init<const G4String&>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double, G4double, G4double, G4double,
	 G4double, G4double, G4double>())
    // ---
    .def("GetZHalfLength",   &G4Trap::GetZHalfLength)
    .def("GetYHalfLength1",  &G4Trap::GetYHalfLength1)
    .def("GetXHalfLength1",  &G4Trap::GetXHalfLength1)
    .def("GetXHalfLength2",  &G4Trap::GetXHalfLength2)
    .def("GetTanAlpha1",     &G4Trap::GetTanAlpha1)
    .def("GetYHalfLength2",  &G4Trap::GetYHalfLength2)
    .def("GetXHalfLength3",  &G4Trap::GetXHalfLength3)
    .def("GetXHalfLength4",  &G4Trap::GetXHalfLength4)
    .def("GetTanAlpha2",     &G4Trap::GetTanAlpha2)
    .def("GetSidePlane",     &G4Trap::GetSidePlane)
    .def("GetSymAxis",       &G4Trap::GetSymAxis)
    .def("SetAllParameters", &G4Trap::SetAllParameters)
    .def("GetCubicVolume",   &G4Trap::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
