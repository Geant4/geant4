//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: pyG4PVPlacement.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
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

