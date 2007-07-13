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
// $Id: pyG4Hype.cc,v 1.2 2007-07-13 04:57:50 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Hype.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Hype.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4Hype {

G4Hype* CreateHype(const G4String& name, 
                   G4double  newInnerRadius,
                   G4double  newOuterRadius,
                   G4double  newInnerStereo,
                   G4double  newOuterStereo,
                   G4double  newHalfLenZ)
{

  return new G4Hype(name, newInnerRadius, newOuterRadius,
                    newInnerStereo, newOuterStereo,
                    newHalfLenZ);
}

}

using namespace pyG4Hype;

// ====================================================================
// module definition
// ====================================================================
void export_G4Hype()
{
  class_<G4Hype, G4Hype*, bases<G4VSolid> >
    ("G4Hype", "hyperbolic solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
                               G4double, G4double>())
    // ---
    .def("GetInnerRadius",   &G4Hype::GetInnerRadius)
    .def("GetOuterRadius",   &G4Hype::GetOuterRadius)
    .def("GetZHalfLength",   &G4Hype::GetZHalfLength)
    .def("GetInnerStereo",   &G4Hype::GetInnerStereo)
    .def("GetOuterStereo",   &G4Hype::GetOuterStereo)
    .def("SetInnerRadius",   &G4Hype::SetInnerRadius)
    .def("SetOuterRadius",   &G4Hype::SetOuterRadius)
    .def("SetZHalfLength",   &G4Hype::SetZHalfLength)
    .def("SetInnerStereo",   &G4Hype::SetInnerStereo)
    .def("SetOuterStereo",   &G4Hype::SetOuterStereo)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateHype", CreateHype, return_value_policy<manage_new_object>());

}

