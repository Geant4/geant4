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
// $Id: pyG4Trd.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Trd.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Trd.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Trd()
{
  class_<G4Trd, G4Trd*, bases<G4VSolid> >
    ("G4Trd", "Trapezoild solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double>())
    // ---
    .def("GetXHalfLength1", &G4Trd::GetXHalfLength1)
    .def("GetXHalfLength2", &G4Trd::GetXHalfLength2)
    .def("GetYHalfLength1", &G4Trd::GetYHalfLength1)
    .def("GetYHalfLength2", &G4Trd::GetYHalfLength2)
    .def("GetZHalfLength",  &G4Trd::GetZHalfLength)
    .def("SetXHalfLength1", &G4Trd::SetXHalfLength1)
    .def("SetXHalfLength2", &G4Trd::SetXHalfLength2)
    .def("SetYHalfLength1", &G4Trd::SetYHalfLength1)
    .def("SetYHalfLength2", &G4Trd::SetYHalfLength2)
    .def("SetZHalfLength",  &G4Trd::SetZHalfLength)
    .def("GetCubicVolume",  &G4Trd::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

