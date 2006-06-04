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
// $Id: pyG4Para.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Para.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Para.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Para()
{
  class_<G4Para, G4Para*, bases<G4VSolid> >
    ("G4Para", "Skewed box sold class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double, G4double>())
    // ---
    .def("GetZHalfLength",    &G4Para::GetZHalfLength)
    .def("GetSymAxis",        &G4Para::GetSymAxis)
    .def("GetYHalfLength",    &G4Para::GetYHalfLength)
    .def("GetXHalfLength",    &G4Para::GetXHalfLength)
    .def("GetTanAlpha",       &G4Para::GetTanAlpha)
    .def("SetXHalfLength",    &G4Para::SetXHalfLength)
    .def("SetYHalfLength",    &G4Para::SetYHalfLength)
    .def("SetZHalfLength",    &G4Para::SetZHalfLength)
    .def("SetAlpha",          &G4Para::SetAlpha)
    .def("SetTanAlpha",       &G4Para::SetTanAlpha)
    .def("SetThetaAndPhi",    &G4Para::SetThetaAndPhi)
    .def("SetAllParameters",  &G4Para::SetAllParameters)
    .def("GetCubicVolume",    &G4Para::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

