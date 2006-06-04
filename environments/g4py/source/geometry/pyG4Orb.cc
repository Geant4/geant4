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
// $Id: pyG4Orb.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Orb.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Orb.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Orb()
{
  class_<G4Orb, G4Orb*, bases<G4VSolid> >
    ("G4Orb", "Orb solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double>())
    // ---
    .def("GetRadius", &G4Orb::GetRadius)
    .def("SetRadius", &G4Orb::SetRadius)
    .def("GetCubicVolume", &G4Orb::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
