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
// $Id: pyG4TwistedTrap.cc,v 1.2 2007-07-13 04:57:50 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4TwistedTrap.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TwistedTrap.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4TwistedTrap {
  
G4TwistedTrap* f1_CreateTwistedTrap(const G4String& name,
                                    G4double  pPhiTwist,
                                    G4double  pDx1, G4double  pDx2,
                                    G4double  pDy, G4double  pDz)
{
  return new G4TwistedTrap(name, pPhiTwist, pDx1, pDx2, pDy, pDz);
}


G4TwistedTrap* f2_CreateTwistedTrap(const G4String& name,
                                    G4double  pPhiTwist,
                                    G4double  pDz, G4double  pTheta,
                                    G4double  pPhi, G4double  pDy1,
                                    G4double  pDx1, G4double  pDx2,
                                    G4double  pDy2, G4double  pDx3,
                                    G4double  pDx4, G4double  pAlph)
{
  return new G4TwistedTrap(name, pPhiTwist, pDz, pTheta, pPhi, 
                           pDy1, pDx1, pDx2, pDy2, pDx3, pDx4, pAlph);
}

}

using namespace pyG4TwistedTrap;

// ====================================================================
// module definition
// ====================================================================
void export_G4TwistedTrap()
{
  class_<G4TwistedTrap, G4TwistedTrap*, bases<G4VSolid> >
    ("G4TwistedTrap", "twisted trapezoid solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
                               G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double,
                               G4double, G4double, G4double,
                               G4double, G4double, G4double,
                               G4double, G4double>())
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateTwistedTap", f1_CreateTwistedTrap,
        return_value_policy<manage_new_object>());
    def("CreateTwistedTap", f2_CreateTwistedTrap,
        return_value_policy<manage_new_object>());

}

