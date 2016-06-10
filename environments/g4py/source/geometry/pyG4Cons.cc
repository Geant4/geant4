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
// $Id: pyG4Cons.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4Cons.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Cons.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4Cons {

G4Cons* CreateCons(const G4String& name, G4double pRmin1, G4double pRmax1,
                   G4double pRmin2, G4double pRmax2, G4double pDz, 
                   G4double pSPhi, G4double pDPhi)
{
  return new G4Cons(name, pRmin1, pRmax1, pRmin2, pRmax2, pDz, pSPhi, pDPhi);
}

}

using namespace pyG4Cons;

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
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateCons", CreateCons, return_value_policy<manage_new_object>());

}

