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
// $Id: pyG4Box.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Box.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Box.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Box()
{
  class_<G4Box, G4Box*, bases<G4VSolid> >
    ("G4Box", "box solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double>())
    // ---
    .def("GetXHalfLength",   &G4Box::GetXHalfLength)
    .def("GetYHalfLength",   &G4Box::GetYHalfLength)
    .def("GetZHalfLength",   &G4Box::GetZHalfLength)
    .def("SetXHalfLength",   &G4Box::SetXHalfLength)
    .def("SetYHalfLength",   &G4Box::SetYHalfLength)
    .def("SetZHalfLength",   &G4Box::SetZHalfLength)
    .def("GetCubicVolume",   &G4Box::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

