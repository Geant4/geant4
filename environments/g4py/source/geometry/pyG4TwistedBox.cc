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
// $Id: pyG4TwistedBox.cc,v 1.2 2007-07-13 04:57:50 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4TwistedBox.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TwistedBox.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4TwistedBox {

G4TwistedBox* CreateTwistedBox(const G4String& name,
                               G4double  pPhiTwist,
                               G4double  pDx, G4double  pDy, G4double  pDz)
{

  return new G4TwistedBox(name, pPhiTwist, pDx, pDy, pDz );

}

}

using namespace pyG4TwistedBox;

// ====================================================================
// module definition
// ====================================================================
void export_G4TwistedBox()
{
  class_<G4TwistedBox, G4TwistedBox*, bases<G4VSolid> >
    ("G4TwistedBox", "twisted box solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    // ---
    .def("GetXHalfLength",   &G4TwistedBox::GetXHalfLength)
    .def("GetYHalfLength",   &G4TwistedBox::GetYHalfLength)
    .def("GetZHalfLength",   &G4TwistedBox::GetZHalfLength)
    .def("GetPhiTwist",      &G4TwistedBox::GetPhiTwist)
    // operators
    .def(self_ns::str(self))
    ;

    // Create solid
    def("CreateTwistedBox", CreateTwistedBox, 
        return_value_policy<manage_new_object>());

}

