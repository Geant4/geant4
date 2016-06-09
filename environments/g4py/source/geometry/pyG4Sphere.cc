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
// $Id: pyG4Sphere.cc,v 1.4 2006/06/29 15:32:27 gunter Exp $
// $Name: geant4-08-01 $
// ====================================================================
//   pyG4Sphere.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Sphere.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Sphere()
{
  class_<G4Sphere, G4Sphere*, bases<G4VSolid> >
    ("G4Sphere", "Sphere solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
         G4double, G4double, G4double>())
    // ---
    .def("GetInsideRadius",     &G4Sphere::GetInsideRadius)
    .def("GetOuterRadius",      &G4Sphere::GetOuterRadius)
    .def("GetStartPhiAngle",    &G4Sphere::GetStartPhiAngle)
    .def("GetDeltaPhiAngle",    &G4Sphere::GetDeltaPhiAngle)
    .def("GetStartThetaAngle",  &G4Sphere::GetStartThetaAngle)
    .def("GetDeltaThetaAngle",  &G4Sphere::GetDeltaThetaAngle)
    .def("SetInsideRadius",     &G4Sphere::SetInsideRadius)
    .def("SetOuterRadius",      &G4Sphere::SetOuterRadius)
    .def("SetStartPhiAngle",    &G4Sphere::SetStartPhiAngle)
    .def("SetDeltaPhiAngle",    &G4Sphere::SetDeltaPhiAngle)
    .def("SetStartThetaAngle",  &G4Sphere::SetStartThetaAngle)
    .def("SetDeltaThetaAngle",  &G4Sphere::SetDeltaThetaAngle)
    .def("GetCubicVolume",      &G4Sphere::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
