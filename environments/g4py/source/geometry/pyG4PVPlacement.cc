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
// $Id: pyG4PVPlacement.cc,v 1.4 2006/06/29 15:32:15 gunter Exp $
// $Name: geant4-09-00 $
// ====================================================================
//   pyG4PVPlacement.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4PVPlacement()
{
  class_<G4PVPlacement, G4PVPlacement*, bases<G4VPhysicalVolume>, 
    boost::noncopyable >
    ("G4PVPlacement", "physical volume placement", no_init)
    // ---
    .def(init<G4RotationMatrix*, const G4ThreeVector&,
	 G4LogicalVolume*, const G4String&,
	 G4LogicalVolume*, G4bool, G4int>())
    .def(init<const G4Transform3D&, G4LogicalVolume*,
	 const G4String&, G4LogicalVolume*, G4bool, G4int>())
    .def(init<G4RotationMatrix*, const G4ThreeVector&,
	 const G4String, G4LogicalVolume*,
	 G4VPhysicalVolume*, G4bool, G4int>())
    .def(init<const G4Transform3D&, const G4String&,
	 G4LogicalVolume*, G4VPhysicalVolume*, G4bool, G4int>())
    ;
}

