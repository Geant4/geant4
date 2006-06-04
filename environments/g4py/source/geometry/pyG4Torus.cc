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
// $Id: pyG4Torus.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Torus.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Torus.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Torus()
{
  class_<G4Torus, G4Torus*, bases<G4VSolid> >
    ("G4Torus", "Torus solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
         G4double, G4double>())
    // ---
    .def("GetRmin", &G4Torus::GetRmin)
    .def("GetRmax", &G4Torus::GetRmax)
    .def("GetRtor", &G4Torus::GetRtor)
    .def("GetSPhi", &G4Torus::GetSPhi)
    .def("GetDPhi", &G4Torus::GetDPhi)
    .def("GetCubicVolume", &G4Torus::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
