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
// $Id: pyG4Para.cc,v 1.5 2007-07-12 10:01:53 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Para.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Para.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4Para {

G4Para* CreatePara(const G4String& name, 
                   G4double pDx, G4double pDy, G4double pDz,
                   G4double pAlpha, G4double pTheta, G4double pPhi)
{
  return new G4Para(name, pDx, pDy, pDz, pAlpha, pTheta, pPhi);
}

}

using namespace pyG4Para;

// ====================================================================
// module definition
// ====================================================================
void export_G4Para()
{
  class_<G4Para, G4Para*, bases<G4VSolid> >
    ("G4Para", "Skewed box sold class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double, G4double>())
    // ---
    .def("GetZHalfLength",    &G4Para::GetZHalfLength)
    .def("GetSymAxis",        &G4Para::GetSymAxis)
    .def("GetYHalfLength",    &G4Para::GetYHalfLength)
    .def("GetXHalfLength",    &G4Para::GetXHalfLength)
    .def("GetTanAlpha",       &G4Para::GetTanAlpha)
    .def("SetXHalfLength",    &G4Para::SetXHalfLength)
    .def("SetYHalfLength",    &G4Para::SetYHalfLength)
    .def("SetZHalfLength",    &G4Para::SetZHalfLength)
    .def("SetAlpha",          &G4Para::SetAlpha)
    .def("SetTanAlpha",       &G4Para::SetTanAlpha)
    .def("SetThetaAndPhi",    &G4Para::SetThetaAndPhi)
    .def("SetAllParameters",  &G4Para::SetAllParameters)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreatePara", CreatePara, return_value_policy<manage_new_object>());
}

