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
// $Id: pyG4PVPlacement.cc,v 1.6 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4PVPlacement.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4PVPlacement {

#if G4VERSION_NUMBER <=711
#elif G4VERSION_NUMBER <= 820
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_CheckOverlaps,
                                       CheckOverlaps, 0, 1);
#elif G4VERSION_NUMBER <=821
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_CheckOverlaps,
                                       CheckOverlaps, 0, 2);
#else
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_CheckOverlaps,
                                       CheckOverlaps, 0, 3);
#endif

}

using namespace pyG4PVPlacement;

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
#if G4VERSION_NUMBER >=800
    .def(init<G4RotationMatrix*, const G4ThreeVector&,
	 G4LogicalVolume*, const G4String&,
	 G4LogicalVolume*, G4bool, G4int, G4bool>())
    .def(init<const G4Transform3D&, G4LogicalVolume*,
	 const G4String&, G4LogicalVolume*, G4bool, G4int, G4bool>())
    .def(init<G4RotationMatrix*, const G4ThreeVector&,
	 const G4String, G4LogicalVolume*,
	 G4VPhysicalVolume*, G4bool, G4int, G4bool>())
    .def(init<const G4Transform3D&, const G4String&,
	 G4LogicalVolume*, G4VPhysicalVolume*, G4bool, G4int, G4bool>())
#endif
    // ---
#if G4VERSION_NUMBER >=800
    .def("CheckOverlaps", &G4PVPlacement::CheckOverlaps, f_CheckOverlaps())
#endif
    ;
}

