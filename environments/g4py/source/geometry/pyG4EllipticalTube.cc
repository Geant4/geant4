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
// $Id: pyG4EllipticalTube.cc,v 1.2 2007-07-13 04:57:50 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4EllipticalTube.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4EllipticalTube.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4EllipticalTube {

G4EllipticalTube* CreateEllipticalTube(const G4String& name,
                                       G4double theDx,
                                       G4double theDy,
                                       G4double theDz)
{
  return new G4EllipticalTube(name, theDx, theDy, theDz);
}

}

using namespace pyG4EllipticalTube;

// ====================================================================
// module definition
// ====================================================================
void export_G4EllipticalTube()
{
  class_<G4EllipticalTube, G4EllipticalTube*, bases<G4VSolid> >
    ("G4EllipticalTube", "elliptical tube solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double>())
    // ---
    .def("GetDx",  &G4EllipticalTube::GetDx)
    .def("GetDy",  &G4EllipticalTube::GetDy)
    .def("GetDz",  &G4EllipticalTube::GetDz)
    .def("SetDx",  &G4EllipticalTube::SetDx)
    .def("SetDy",  &G4EllipticalTube::SetDy)
    .def("SetDz",  &G4EllipticalTube::SetDz)

    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateEllipticalTube", CreateEllipticalTube,
        return_value_policy<manage_new_object>());

}

