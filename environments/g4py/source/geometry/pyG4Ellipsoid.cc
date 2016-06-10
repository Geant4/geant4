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
// $Id: pyG4Ellipsoid.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4Ellipsoid.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Ellipsoid.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Ellipsoid {

G4Ellipsoid* CreateEllipsoid(const G4String& name,
                             G4double  pxSemiAxis,
                             G4double  pySemiAxis,
                             G4double  pzSemiAxis,
                             G4double  pzBottomCut=0,
                             G4double  pzTopCut=0)
{
  return new G4Ellipsoid(name, pxSemiAxis, pySemiAxis, pzSemiAxis,
                         pzBottomCut, pzTopCut);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(f_CreateEllipsoid, CreateEllipsoid, 4,6)

}

using namespace pyG4Ellipsoid;

// ====================================================================
// module definition
// ====================================================================
void export_G4Ellipsoid()
{
  class_<G4Ellipsoid, G4Ellipsoid*, bases<G4VSolid> >
    ("G4Ellipsoid", "ellipsoid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double,
                               G4double>())
    // ---
    .def("GetSemiAxisMax", &G4Ellipsoid::GetSemiAxisMax)
    .def("GetZBottomCut",  &G4Ellipsoid::GetZBottomCut)
    .def("GetZTopCut",     &G4Ellipsoid::GetZTopCut)
    .def("SetSemiAxis",  &G4Ellipsoid::SetSemiAxis)
    .def("SetZCuts",  &G4Ellipsoid::SetZCuts)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateEllipsoid", CreateEllipsoid,
        f_CreateEllipsoid()[return_value_policy<manage_new_object>()]);
}

