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
// $Id: pyG4TwistedTrd.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4TwistedTrap.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TwistedTrd.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4TwistedTrd {

G4TwistedTrd* CreateTwistedTrd(const G4String& name,
                               G4double  pDx1, G4double  pDx2,
                               G4double  pDy1, G4double  pDy2,
                               G4double  pDz, G4double  pPhiTwist)
{
  return new G4TwistedTrd(name, pDx1, pDx2, pDy1, pDy2, pDz, pPhiTwist);
}

}

using namespace pyG4TwistedTrd;

// ====================================================================
// module definition
// ====================================================================
void export_G4TwistedTrd()
{
  class_<G4TwistedTrd, G4TwistedTrd*, bases<G4VSolid> >
    ("G4TwistedTrd", "twisted trapezoid solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
                    G4double, G4double, G4double>())
    // ---
    .def("GetX1HalfLength",    &G4TwistedTrd::GetX1HalfLength)
    .def("GetX2HalfLength",    &G4TwistedTrd::GetX2HalfLength)
    .def("GetY1HalfLength",    &G4TwistedTrd::GetY1HalfLength)
    .def("GetY2HalfLength",    &G4TwistedTrd::GetY2HalfLength)
    .def("GetZHalfLength",     &G4TwistedTrd::GetZHalfLength)
    .def("GetPhiTwist",        &G4TwistedTrd::GetPhiTwist)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateTwistedTrd", CreateTwistedTrd,
        return_value_policy<manage_new_object>());

}

