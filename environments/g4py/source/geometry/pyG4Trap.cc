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
// $Id: pyG4Trap.cc,v 1.6 2007-07-13 04:57:50 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Trap.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Trap.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4Trap {

G4Trap* f1_CreateTrap(const G4String& name)
{
  return new G4Trap(name);
}


G4Trap* f2_CreateTrap(const G4String& name, G4double pDz,
                      G4double pTheta, G4double pPhi,
                      G4double pDy1, G4double pDx1, G4double pDx2,
                      G4double pAlp1,
                      G4double pDy2, G4double pDx3, G4double pDx4,
                      G4double pAlp2)
{
  return new G4Trap(name, pDz, pTheta, pPhi,
                    pDy1, pDx1, pDx2, pAlp1,
                    pDy2, pDx3, pDx4, pAlp2);

}


G4Trap* f3_CreateTrap(const G4String& name,
                      const std::vector<G4ThreeVector>& pt)
{
  G4ThreeVector ptlist[8];
  for (G4int i=0; i<8; i++) {
    ptlist[i]= pt[i];
  }

  return new G4Trap(name, ptlist);
}


G4Trap* f4_CreateTrap(const G4String& name, G4double pZ, 
                      G4double pY, G4double pX, G4double pLTX)
{
  return new G4Trap(name, pZ, pY, pX, pLTX);
}


G4Trap* f5_CreateTrap(const G4String& name, G4double pDx1,  G4double pDx2,
                      G4double pDy1,  G4double pDy2, G4double pDz)
{
  return new G4Trap(name, pDx1, pDx2, pDy1, pDy2, pDz);
}


G4Trap* f6_CreateTrap(const G4String& name, G4double pDx, G4double pDy,
                      G4double pDz, 
                      G4double pAlpha, G4double pTheta, G4double pPhi )
{
  return new G4Trap(name, pDx, pDy, pDz, pAlpha, pTheta, pPhi);
}

}

using namespace pyG4Trap;

// ====================================================================
// module definition
// ====================================================================
void export_G4Trap()
{
  class_<G4Trap, G4Trap*, bases<G4VSolid> >
    ("G4Trap", "Generic trapezoild soild class", no_init)
    // constructors
    .def(init<const G4String&>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double, G4double, G4double, G4double,
	 G4double, G4double, G4double>())
    // ---
    .def("GetZHalfLength",   &G4Trap::GetZHalfLength)
    .def("GetYHalfLength1",  &G4Trap::GetYHalfLength1)
    .def("GetXHalfLength1",  &G4Trap::GetXHalfLength1)
    .def("GetXHalfLength2",  &G4Trap::GetXHalfLength2)
    .def("GetTanAlpha1",     &G4Trap::GetTanAlpha1)
    .def("GetYHalfLength2",  &G4Trap::GetYHalfLength2)
    .def("GetXHalfLength3",  &G4Trap::GetXHalfLength3)
    .def("GetXHalfLength4",  &G4Trap::GetXHalfLength4)
    .def("GetTanAlpha2",     &G4Trap::GetTanAlpha2)
    .def("GetSidePlane",     &G4Trap::GetSidePlane)
    .def("GetSymAxis",       &G4Trap::GetSymAxis)
    .def("SetAllParameters", &G4Trap::SetAllParameters)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateTrap", f1_CreateTrap, return_value_policy<manage_new_object>());
    def("CreateTrap", f2_CreateTrap, return_value_policy<manage_new_object>());
    def("CreateTrap", f3_CreateTrap, return_value_policy<manage_new_object>());
    def("CreateTrap", f4_CreateTrap, return_value_policy<manage_new_object>());
    def("CreateTrap", f5_CreateTrap, return_value_policy<manage_new_object>());
    def("CreateTrap", f6_CreateTrap, return_value_policy<manage_new_object>());

}

