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
// $Id: pyG4EllipticalCone.cc 81291 2014-05-26 09:31:19Z gcosmo $
// ====================================================================
//   pyG4EllipticalCone.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4EllipticalCone.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4EllipticalCone {

G4EllipticalCone* CreateEllipticalCone(const G4String& name, 
                                       G4double  pxSemiAxis,
                                       G4double  pySemiAxis,
                                       G4double  zMax,
                                       G4double  pzTopCut)
{
  return new G4EllipticalCone(name, pxSemiAxis,pySemiAxis, zMax, pzTopCut);
}

}

using namespace pyG4EllipticalCone;

// ====================================================================
// module definition
// ====================================================================
void export_G4EllipticalCone()
{
  class_<G4EllipticalCone, G4EllipticalCone*, bases<G4VSolid> >
    ("G4EllipticalCone", "elliptical cone solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    // ---
    .def("GetSimiAxisMax",  &G4EllipticalCone::GetSemiAxisMax)
    .def("GetZTopCut",      &G4EllipticalCone::GetZTopCut)
    .def("SetSemiAxis",     &G4EllipticalCone::SetSemiAxis)
    .def("SetZCut",         &G4EllipticalCone::SetZCut)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateEllipticalCone", CreateEllipticalCone,
        return_value_policy<manage_new_object>());

}

